Shader "Hidden/FXAA" {
	Properties {
		_MainTex ("Texture", 2D) = "white" {}
	}

	CGINCLUDE
		#include "UnityCG.cginc"

		sampler2D _MainTex;
		float4 _MainTex_TexelSize;

		float _ContrastThreshold, _RelativeThreshold;
		float _SubpixelBlending;

		struct VertexData {
			float4 vertex : POSITION;
			float2 uv : TEXCOORD0;
		};

		struct Interpolators {
			float4 pos : SV_POSITION;
			float2 uv : TEXCOORD0;
		};

		Interpolators VertexProgram (VertexData v) {
			Interpolators i;
			i.pos = UnityObjectToClipPos(v.vertex);
			i.uv = v.uv;
			return i;
		}

		float4 Sample (float2 uv) {
			return tex2Dlod(_MainTex, float4(uv, 0, 0));
		}

		float SampleLuminance (float2 uv) {
			#if defined(LUMINANCE_GREEN)
				return Sample(uv).g;
			#else
				return Sample(uv).a;
			#endif
		}

		float SampleLuminance (float2 uv, float uOffset, float vOffset) {
			uv += _MainTex_TexelSize * float2(uOffset, vOffset);
			return SampleLuminance(uv);
		}

		struct LuminanceData {
			float m, n, e, s, w;
			float ne, nw, se, sw;
			float highest, lowest, contrast;
		};
        //一次性采样所有需要用到的9个亮度值，并计算contrast
		LuminanceData SampleLuminanceNeighborhood (float2 uv) {
			LuminanceData l;
			l.m = SampleLuminance(uv);
			l.n = SampleLuminance(uv,  0,  1);
			l.e = SampleLuminance(uv,  1,  0);
			l.s = SampleLuminance(uv,  0, -1);
			l.w = SampleLuminance(uv, -1,  0);

			l.ne = SampleLuminance(uv,  1,  1);
			l.nw = SampleLuminance(uv, -1,  1);
			l.se = SampleLuminance(uv,  1, -1);
			l.sw = SampleLuminance(uv, -1, -1);

			l.highest = max(max(max(max(l.n, l.e), l.s), l.w), l.m);//五个值中亮度最大的
			l.lowest = min(min(min(min(l.n, l.e), l.s), l.w), l.m);//五个值中亮度最小的
			l.contrast = l.highest - l.lowest;
			return l;
		}

		bool ShouldSkipPixel (LuminanceData l) {
			float threshold =
				max(_ContrastThreshold, _RelativeThreshold * l.highest);
			return l.contrast < threshold;
		}
        //确定像素的BlendFactor
		float DeterminePixelBlendFactor (LuminanceData l) {
			float filter = 2 * (l.n + l.e + l.s + l.w);
			filter += l.ne + l.nw + l.se + l.sw;
			filter *= 1.0 / 12;
			filter = abs(filter - l.m);
			filter = saturate(filter / l.contrast);

			float blendFactor = smoothstep(0, 1, filter);
			return blendFactor * blendFactor * _SubpixelBlending;
		}

		struct EdgeData {
			bool isHorizontal;
			float pixelStep;
			float oppositeLuminance, gradient;
		};
        //确定边的走向，其实是确定的九宫格中9个采样值所组成的色块的边沿的走向。
		EdgeData DetermineEdge (LuminanceData l) {
			EdgeData e;
			float horizontal =
				abs(l.n + l.s - 2 * l.m) * 2 +
				abs(l.ne + l.se - 2 * l.e) +
				abs(l.nw + l.sw - 2 * l.w);
			float vertical =
				abs(l.e + l.w - 2 * l.m) * 2 +
				abs(l.ne + l.nw - 2 * l.n) +
				abs(l.se + l.sw - 2 * l.s);
			e.isHorizontal = horizontal >= vertical;//通过计算水平和垂直方向不同权重来判断走向。

			float pLuminance = e.isHorizontal ? l.n : l.e;//根据走向判断反方向的亮度值需要取哪个值。
			float nLuminance = e.isHorizontal ? l.s : l.w;
			float pGradient = abs(pLuminance - l.m);//计算正方向的亮度梯度
			float nGradient = abs(nLuminance - l.m);//计算反方向的亮度梯度

			e.pixelStep =
				e.isHorizontal ? _MainTex_TexelSize.y : _MainTex_TexelSize.x;//根据方向确定像素步长
			
			if (pGradient < nGradient) {//根据正向和反向的亮度梯度，判断最终blend的方向。
				e.pixelStep = -e.pixelStep;
				e.oppositeLuminance = nLuminance;
				e.gradient = nGradient;
			}
			else {
				e.oppositeLuminance = pLuminance;
				e.gradient = pGradient;
			}

			return e;
		}

		#if defined(LOW_QUALITY)//低质量的情况下沿边沿方向查看4的采样点来确定端点；
			#define EDGE_STEP_COUNT 4
			#define EDGE_STEPS 1, 1.5, 2, 4//采样点的uv增量每次都不同
			#define EDGE_GUESS 12
		#else//高质量的情况下沿边沿方向查看10的采样点来确定端点；
			#define EDGE_STEP_COUNT 10
			#define EDGE_STEPS 1, 1.5, 2, 2, 2, 2, 2, 2, 2, 4//采样点的uv增量每次都不同
			#define EDGE_GUESS 8
		#endif

		static const float edgeSteps[EDGE_STEP_COUNT] = { EDGE_STEPS };//构造一个数组

        //计算边沿的BlenderFactor
		float DetermineEdgeBlendFactor (LuminanceData l, EdgeData e, float2 uv) {
			float2 uvEdge = uv;
			float2 edgeStep;
			if (e.isHorizontal) {
				uvEdge.y += e.pixelStep * 0.5;
				edgeStep = float2(_MainTex_TexelSize.x, 0);
			}
			else {
				uvEdge.x += e.pixelStep * 0.5;
				edgeStep = float2(0, _MainTex_TexelSize.y);
			}

			float edgeLuminance = (l.m + e.oppositeLuminance) * 0.5;
			float gradientThreshold = e.gradient * 0.25;

			float2 puv = uvEdge + edgeStep * edgeSteps[0];
			float pLuminanceDelta = SampleLuminance(puv) - edgeLuminance;
			bool pAtEnd = abs(pLuminanceDelta) >= gradientThreshold;

			UNITY_UNROLL//告诉编译器编译成非循环的形式
			for (int i = 1; i < EDGE_STEP_COUNT && !pAtEnd; i++) {
				puv += edgeStep * edgeSteps[i];//顺边沿方向运动采样坐标
				pLuminanceDelta = SampleLuminance(puv) - edgeLuminance;//计算采样的亮度值与起始点的亮度值的差
				pAtEnd = abs(pLuminanceDelta) >= gradientThreshold;//判断是否大于起始点的亮度梯度
			}
			if (!pAtEnd) {//没有找到结尾点的情况下，猜测一个。
				puv += edgeStep * EDGE_GUESS;
			}

			float2 nuv = uvEdge - edgeStep * edgeSteps[0];
			float nLuminanceDelta = SampleLuminance(nuv) - edgeLuminance;
			bool nAtEnd = abs(nLuminanceDelta) >= gradientThreshold;

			UNITY_UNROLL
			for (int i = 1; i < EDGE_STEP_COUNT && !nAtEnd; i++) {
				nuv -= edgeStep * edgeSteps[i];
				nLuminanceDelta = SampleLuminance(nuv) - edgeLuminance;
				nAtEnd = abs(nLuminanceDelta) >= gradientThreshold;
			}
			if (!nAtEnd) {
				nuv -= edgeStep * EDGE_GUESS;
			}

			float pDistance, nDistance;
			if (e.isHorizontal) {
				pDistance = puv.x - uv.x;
				nDistance = uv.x - nuv.x;
			}
			else {
				pDistance = puv.y - uv.y;
				nDistance = uv.y - nuv.y;
			}

			float shortestDistance;
			bool deltaSign;
			if (pDistance <= nDistance) {
				shortestDistance = pDistance;
				deltaSign = pLuminanceDelta >= 0;
			}
			else {
				shortestDistance = nDistance;
				deltaSign = nLuminanceDelta >= 0;
			}

			if (deltaSign == (l.m - edgeLuminance >= 0)) {
				return 0;
			}
			return 0.5 - shortestDistance / (pDistance + nDistance);
		}

		float4 ApplyFXAA (float2 uv) {
			LuminanceData l = SampleLuminanceNeighborhood(uv);//获取亮度信息
			if (ShouldSkipPixel(l)) {//如果Contrast小于标准，直接返回像素值，不需要抗锯齿。
				return Sample(uv);
			}

			float pixelBlend = DeterminePixelBlendFactor(l);//计算像素的BlenderFactor
			EdgeData e = DetermineEdge(l);//判断边沿的走向
			float edgeBlend = DetermineEdgeBlendFactor(l, e, uv);//计算边沿的BlenderFactor
			float finalBlend = max(pixelBlend, edgeBlend);//取最大值作为最终的BlenderFactor

			if (e.isHorizontal) {//根据边沿的走向和BlenderFactor确定最终的偏移uv
				uv.y += e.pixelStep * finalBlend;
			}
			else {
				uv.x += e.pixelStep * finalBlend;
			}
			return float4(Sample(uv).rgb, l.m);//按照偏移之后的uv值采样后，返回这个采样值。同时返回的中间点的亮度值以备他用。
		}
	ENDCG

	SubShader {
		Cull Off
		ZTest Always
		ZWrite Off

		Pass { // 0 luminancePass
			CGPROGRAM
				#pragma vertex VertexProgram
				#pragma fragment FragmentProgram

				#pragma multi_compile _ GAMMA_BLENDING

				float4 FragmentProgram (Interpolators i) : SV_Target {
					float4 sample = tex2D(_MainTex, i.uv);
					sample.rgb = saturate(sample.rgb);//HDR转换到LDR
					sample.a = LinearRgbToLuminance(sample.rgb);//通过颜色计算亮度的函数
					#if defined(GAMMA_BLENDING)
						sample.rgb = LinearToGammaSpace(sample.rgb);//将颜色从线性空间转到gamma空间，用于FXAA的计算。因为在线性空间下计算出得结果虽然不准确但是视觉效果好。
					#endif
					return sample;
				}
			ENDCG
		}

		Pass { // 1 fxaaPass
			CGPROGRAM
				#pragma vertex VertexProgram
				#pragma fragment FragmentProgram

				#pragma multi_compile _ LUMINANCE_GREEN
				#pragma multi_compile _ LOW_QUALITY//在DetermineEdgeBlendFactor时向两侧排查的距离只有原来的一半。
				#pragma multi_compile _ GAMMA_BLENDING

				float4 FragmentProgram (Interpolators i) : SV_Target {
					float4 sample = ApplyFXAA(i.uv);
					#if defined(GAMMA_BLENDING)
						sample.rgb = GammaToLinearSpace(sample.rgb);//将颜色从Gamma空间转到线性空间
					#endif
					return sample;
				}
			ENDCG
		}
	}
}