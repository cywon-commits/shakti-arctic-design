# 선행 연구 점검: 오목한 도미노 영역 (2026-11-25)

질문: 십자 시료의 문턱(c = n/3, c = n/√2)과 class 10 곡선이 이미 알려져 있는가.
검색은 웹 검색과 arXiv 초록·본문 일부를 읽은 수준이다. 본문 전체를 읽지 않은 논문은 그렇게 표시했다.

## 요약

1. **모서리 하나를 자른 Aztec diamond는 잘 연구되어 있다.** 문턱, 곡선, 요동이 모두 나와 있다.
2. **네 모서리를 모두 자른 십자 모양은 찾지 못했다.** 검색한 범위에서는 c/n의 닫힌 식도, 타일링이 사라지는 문턱 c = n/3도 없었다.
3. **곡선이 대수 곡선이라는 것은 일반 정리다.** 따라서 "유리 곡선으로 풀린다"는 것 자체는 새롭지 않다. 새로운 것은 구체적인 식과 Shakti로의 번역이다.
4. **stem-island 규칙, 가장자리 사전, 모서리 절단 규칙은 Shakti 고유의 내용이라 대응하는 선행 연구가 없다.**

## 논문별

| 논문 | 영역 | 결과 | 우리 결과와의 관계 |
|---|---|---|---|
| Colomo, Pronko, arXiv:1306.6207 | 한 모서리에서 정사각형을 떼어낸 Aztec diamond | 잘린 부분이 원래의 arctic 타원에 닿을 때 자유에너지의 3차 상전이. v_c = s/r = (1−√α)/(1+√α) | 십자의 c = n/√2와 **같은 조건**. α = 1/2에서 s/(r+s) = (1−1/√2)/2이고, 십자에서는 (n−c)/(2n)이 이 값이다. 일치 확인. |
| Colomo, Pronko, Sportiello, arXiv:1807.07549 (J. Stat. Phys. 2018) | L자 영역 (free-fermion 6-vertex, DWBC) | Tangent method로 변형된 arctic 곡선의 매개변수 식. 경계와의 접점 6개, 첨점 2개 | 오목 모서리 하나의 경우. 우리 곡선의 "오목 모서리마다 생기는 주머니"에 해당한다. 4절 뒷부분과 결론은 읽지 못했다. |
| Colomo, Pronko, arXiv:2407.07849 (PRE 110, 054140) | 삼각형을 떼어낸 오각형 영역 | 다섯째 변이 타원에 접할 때 3차 상전이. 곡선 자체는 구하지 않음 | 45° 절단 하나에 해당. 우리 모서리 절단 규칙(볼록 절단은 원을 바꾸지 않음)과 상보적이다. |
| Charlier, Claeys, arXiv:2512.20388 | L자 Aztec 영역 | 타일링 수의 정밀한 점근식, 네 영역, 임계점에서 Tracy–Widom | 문턱 근처의 유한 크기 보정이 필요하면 여기서 가져올 수 있다. |
| Ferrari, Vető, arXiv:1909.10840 (Ann. Appl. Probab.) | 위쪽을 수평선으로 자른 Aztec diamond | 절단선이 원에 접할 때 hard-edge tacnode 과정 | c = n/√2 바로 위에서의 경계 요동. |
| Astala, Duse, Prause, Zhong, arXiv:2004.02599 | 다각형 영역 (일반 dimer 모형) | 얼어붙은 경계는 대수 곡선, 일반점에서 Pokrovsky–Talapov | 십자 곡선이 대수적이라는 것의 일반 근거. 오목 다각형이 포함되는지는 초록만으로 확인하지 못했다. |
| Kenyon, Prause, arXiv:2006.01219; Kuchumov, arXiv:2603.21255 | tangent plane method; 구멍 뚫린 Aztec diamond | 타원함수 매개변수 식 | 우리가 쓴 방법(접선족)의 표준적인 형태. 구멍 있는 시료로 넓힐 때의 출발점. |

읽지 못하고 제목만 확인한 것: "Boundary dents, the arctic circle and the arctic ellipse" (arXiv:2308.06863), Bufetov–Knizel의 rectangular Aztec diamond, Di Francesco–Guitter의 tangent method 논문들.

## 연구 방향에 주는 뜻

- **주장의 무게 중심을 옮긴다.** "오목 영역의 arctic 곡선을 풀었다"는 수학적 신규성으로 내세우기 어렵다. 내세울 것은 (i) Shakti 시료의 절단·종단이 어떤 도미노 영역을 만드는지의 사전, (ii) 그 사전으로 액체 영역의 위치·크기·모양을 설계한다는 것, (iii) 십자에 대한 닫힌 식이다.
- **c = n/√2는 3차 상전이다.** 지금까지는 "곡선이 원에서 벗어나기 시작하는 곳"으로만 적었다. 선행 연구에 따르면 바닥상태 엔트로피의 3차 도함수가 뛴다. Shakti에서는 팔 폭을 바꿀 때 잔류 엔트로피에 나타나는 3차 특이성이 된다. 인용하고, 네 모서리에서도 3차인지 확인할 가치가 있다.
- **오목 모서리 조건을 검증했다.** 같은 접선족 풀이를 L자 영역에 적용하면 Colomo–Pronko–Sportiello의 식 (4.15)가 재현된다 (`analytic-results.md` 2.6절).
- **c = n/3은 높이 함수의 Lipschitz 조건에서 바로 나오므로** 알려져 있지 않더라도 정리로 내세울 내용은 아니다.

## 참고

- https://arxiv.org/abs/1306.6207
- https://arxiv.org/abs/1807.07549
- https://arxiv.org/abs/2407.07849
- https://arxiv.org/abs/2512.20388
- https://arxiv.org/abs/1909.10840
- https://arxiv.org/abs/2004.02599
- https://arxiv.org/abs/2006.01219
- https://arxiv.org/abs/2603.21255
