# Mini S3 - S3-inspired Distributed Storage Simulator in C

AWS S3 같은 대규모 오브젝트 스토리지에서 자주 등장하는 개념을 작은 C 프로그램으로 실험하는 교육용 시뮬레이터.

이 프로젝트는 다음 3가지를 직접 구현하고 검증한다.

1. 데이터 분산 저장 (data striping)
2. Power of Two Choices 기반 단순 로드 밸런싱
3. Reed-Solomon 이레이저 코딩 (Cauchy matrix + GF(256))

이 레포는 실제 AWS S3 구현을 재현한 프로젝트가 아니다. 단일 프로세스, 단일 머신 환경에서 분산 저장의 핵심 아이디어를 이해하고 실험하기 위한 toy model이다.

## What This Project Implements

### 1. Data Striping
데이터를 여러 조각으로 나누고 서로 다른 디스크에 저장하는 개념을 데모로 보여준다. 실제 분산 시스템의 네트워크/디스크 병렬성까지 재현하지는 않지만, "조각 단위 병렬 처리"라는 아이디어를 설명하기에는 충분하다.

### 2. Power of Two Choices
매번 임의의 디스크 2개를 고른 뒤 더 한가한 쪽에 배치한다. 구현은 매우 단순하지만, 순수 랜덤 배치에 비해 최대 부하와 표준편차가 얼마나 줄어드는지 통계적으로 확인할 수 있다.

### 3. Erasure Coding
`k=5` 데이터 조각과 `m=4` 패리티 조각을 사용하는 Reed-Solomon 코드를 구현했다. 총 9개 조각 중 아무 5개만 남아 있으면 원본 데이터를 복구할 수 있다. 이 설정은 이 프로젝트의 실험용 파라미터이며, 특정 상용 시스템의 내부 설정을 주장하지 않는다.

## What This Project Is Not

- 실제 S3 내부 구현 공개판이 아니다.
- 멀티노드 스토리지 서버, 네트워크 스택, 디스크 스케줄러, 복제 파이프라인을 포함하지 않는다.
- 실제 운영 환경의 상관 고장, rack/AZ 배치, 재해, 소프트웨어 버그, 운영 실수까지 모델링하지 않는다.
- 처리량 수치는 이 코드와 현재 머신에서 측정한 예시일 뿐이며, 환경이 바뀌면 달라진다.

## Build & Run

```bash
make        # 빌드
make run    # 3가지 개념 데모 실행
make test   # 검증/실험 실행
make clean  # 정리
```

## File Layout

```text
gf256.h / gf256.c       - GF(256) finite field arithmetic
erasure.h / erasure.c   - Cauchy matrix based Reed-Solomon erasure coding
main.c                  - striping / P2C / erasure coding demo
test.c                  - validation tests and simulation experiments
Makefile                - build rules
```

## Validation And Experiments

`make test`는 기능 검증과 toy simulation 결과를 함께 출력한다. 난수 기반 실험은 `srand(42)`로 고정되어 있어 확률 실험 흐름은 재현 가능하다. 다만 wall-clock 기반 처리량 수치는 머신, 부하, 컴파일러 옵션에 따라 달라질 수 있다.

### 1. GF(256) Field Checks

GF(256) 구현이 기본 필드 성질을 만족하는지 검사한다.

| Check | Result |
|------|------|
| Multiplicative identity (`a * 1 = a`) | 256/256 |
| Multiplicative inverse (`a * a^-1 = 1`) | 255/255 |
| Commutativity (`a * b = b * a`) | 65536/65536 |
| Distributivity (`a * (b + c) = ab + ac`) | 100000/100000 |

### 2. Exhaustive Erasure-Coding Recovery

`RS(5,4)` 설정에서 가능한 모든 1-4개 조각 고장 조합을 전수 검사한다.

| Failed shards | Patterns | Result |
|-------------|---------|------|
| 1 | C(9,1) = 9 | 9/9 |
| 2 | C(9,2) = 36 | 36/36 |
| 3 | C(9,3) = 84 | 84/84 |
| 4 | C(9,4) = 126 | 126/126 |
| Total | 255 | 255/255 PASS |

### 3. Random Recovery Tests

서로 다른 데이터 크기에 대해 랜덤 1-4개 고장을 주고 복구를 반복 검증한다.

| Data size | Trials | Result |
|----------|------|------|
| 100 B | 100 | 100/100 |
| 1 KB | 100 | 100/100 |
| 64 KB | 50 | 50/50 |
| 1 MB | 20 | 20/20 |

### 4. Power of Two Choices Statistics

디스크당 평균 100개 아이템이 들어가도록 설정하고 200회 반복 측정한다.

| Disks | Strategy | Avg max | Avg min | Stddev | Max/mean |
|------|------|-------:|-------:|------:|--------:|
| 16 | rand | 118.2 | 82.7 | 9.5 | 1.1816 |
| 16 | P2C | 101.2 | 98.2 | 0.9 | 1.0122 |
| 64 | rand | 124.2 | 77.2 | 9.9 | 1.2419 |
| 64 | P2C | 101.8 | 97.3 | 0.9 | 1.0175 |
| 256 | rand | 129.8 | 73.0 | 10.0 | 1.2983 |
| 256 | P2C | 102.0 | 96.3 | 1.0 | 1.0202 |
| 512 | rand | 131.9 | 71.0 | 10.0 | 1.3185 |
| 512 | P2C | 102.0 | 96.0 | 1.0 | 1.0203 |

핵심 해석: 이 toy model에서는 디스크 수가 커질수록 순수 랜덤 배치의 최대 부하 편차가 커지지만, P2C는 매우 좁은 범위에 부하를 유지한다.

### 5. Encoding / Decoding Throughput

다음 표는 한 머신에서 측정한 예시 결과다. 절대 성능 수치보다 인코딩/디코딩 경향을 보는 용도로 해석하는 편이 맞다.

Apple M1 Pro 기준. "before"는 최적화 전(byte-at-a-time 루프 + 스칼라), "after"는 최적화 후(shard-major 루프 + NEON SIMD). 총 개선은 루프 구조 변경(~2x)과 SIMD(~17x) 두 가지의 합산이다. SIMD 단독 기여는 [11. SIMD Optimization](#11-simd-optimization) 참조.

| Data size | Enc. before | Enc. after | Enc. speedup | Dec. before | Dec. after | Dec. speedup |
|----------|------------:|-----------:|-------------:|------------:|-----------:|-------------:|
| 1 KB | ~195 | ~1,794 | 9.2x | ~242 | ~1,263 | 5.2x |
| 4 KB | ~199 | ~4,490 | 22.6x | ~265 | ~3,762 | 14.2x |
| 16 KB | ~199 | ~5,914 | 29.7x | ~264 | ~6,191 | 23.5x |
| 64 KB | ~193 | ~6,372 | 33.0x | ~266 | ~7,298 | 27.4x |
| 256 KB | ~196 | ~6,428 | 32.8x | ~267 | ~7,606 | 28.5x |
| 1 MB | ~190 | ~6,874 | 36.2x | ~264 | ~8,571 | 32.5x |

단위: MB/s.

### 6. Tail-Latency Simulation

정상 디스크는 1-5ms, 느린 디스크는 5% 확률로 50-200ms가 걸린다고 가정한다. 비교 대상은 다음 두 가지다.

- naive read: 처음 5개 데이터 조각만 읽고 가장 느린 응답을 기다림
- hedged read: 9개 조각을 모두 요청하고 가장 빠른 5개만 모이면 완료

| Percentile | Naive | Hedged | Improvement |
|-----------|------:|------:|------------:|
| P50 | 5.0 ms | 3.0 ms | 40.0% |
| P90 | 137.0 ms | 4.0 ms | 97.1% |
| P95 | 169.0 ms | 4.0 ms | 97.6% |
| P99 | 195.0 ms | 5.0 ms | 97.4% |
| P99.9 | 200.0 ms | 5.0 ms | 97.5% |

이 결과는 "여분 조각이 있으면 tail latency를 우회할 수 있다"는 개념을 설명한다. 실제 서비스의 네트워크, 큐잉, 백엔드 분산 구조까지 반영한 결과는 아니다.

### 7. Storage-Efficiency Comparison

| Scheme | Storage overhead | Tolerated disk failures | Minimum shards |
|-------|-----------------:|------------------------:|---------------:|
| No redundancy | 1.0x | 0 | 1 |
| 2-way replication | 2.0x | 1 | 2 |
| 3-way replication | 3.0x | 2 | 3 |
| RS(5,4) | 1.8x | 4 | 9 |
| RS(10,4) | 1.4x | 4 | 14 |
| RS(16,4) | 1.25x | 4 | 20 |

이 비교는 purely combinational storage trade-off를 보여준다. 실제 시스템 설계는 repair bandwidth, placement policy, correlated failure까지 같이 봐야 한다.

### 8. Durability Analysis

`test.c`는 단순화된 연속시간 고장/복구 모델을 사용해 MTTDL과 연간 손실 확률을 계산한다.

가정:

- 디스크 고장은 독립적이다.
- AFR은 고정이다.
- rebuild 시간은 일정하다.
- 데이터 손실 조건은 "허용 가능한 고장 수 초과"뿐이다.

이 가정 아래에서 `AFR=2%`, `rebuild=24h`일 때 다음 값을 계산한다.

| Scheme | MTTDL (years) | Annual loss probability | Nines |
|-------|--------------:|------------------------:|------:|
| Single disk | 4.95e+01 | 2.00e-02 | 1.7 |
| 2-way replication | 4.47e+05 | 2.24e-06 | 5.7 |
| 3-way replication | 5.39e+09 | 1.86e-10 | 9.7 |
| RS(5,4) | 8.26e+15 | 1.21e-16 | 15.9 |
| RS(10,4) | 5.31e+14 | 1.88e-15 | 14.7 |

중요: 이 숫자는 이 toy model의 결과다. 실제 클라우드 스토리지의 내구성은 디스크 고장 독립성만으로 설명되지 않는다.

Monte Carlo 검증도 함께 수행해 해석식과 시뮬레이션이 대략 일치하는지 확인한다.

| Scheme | Analytic | Monte Carlo | Relative error |
|-------|---------:|------------:|---------------:|
| Single disk | 20.0000% | 20.0152% | 0.1% |
| 2-way replication | 0.0272% | 0.0245% | 10.0% |
| RS(2,1) | 0.0816% | 0.0804% | 1.4% |

### 9. CRC32 Silent-Corruption Check

고장이 아니라 "틀린 데이터를 정상처럼 반환하는" 상황을 단순하게 시뮬레이션한다. 각 조각에 CRC32를 붙여 무결성을 확인하고, 손상된 조각은 unavailable로 처리한 뒤 이레이저 코딩으로 복구한다.

| Check | Result |
|------|------|
| CRC32 detection rate | 999/1000 (99.9%) |
| Recovery after detection | 1000/1000 (100%) |

### 10. Rebuild Vulnerability Window

디스크 하나가 이미 고장 난 상태에서 rebuild 중 추가 고장이 발생할 확률을 계산한다.

| Rebuild time | 3-way replication | RS(5,4) | Ratio |
|-------------|------------------:|---------:|------:|
| 1h | 5.32e-12 | 1.98e-21 | 2,685,933,369 |
| 6h | 1.91e-10 | 2.57e-18 | 74,612,874 |
| 24h | 3.06e-09 | 6.57e-16 | 4,664,118 |
| 72h | 2.76e-08 | 5.32e-14 | 518,476 |
| 168h | 1.50e-07 | 1.57e-12 | 95,319 |
| 720h | 2.75e-06 | 5.28e-10 | 5,217 |

이 표도 같은 단순화 가정에 의존한다. 의미 있는 메시지는 "rebuild window가 길수록 redundancy scheme 차이가 커진다"는 점이다.

### 11. SIMD Optimization

GF(256) 곱셈-누적(`dst[i] ^= c * src[i]`)은 이레이저 코딩의 핫 패스다. 스칼라 구현은 바이트마다 exp/log 테이블 룩업을 하지만, SIMD는 split-table 기법을 사용한다.

**Split-table 기법:**
각 바이트 `b`를 high/low nibble로 분리한 뒤, 16-entry 테이블 2개로 병렬 lookup:

```
c * b = lo_tbl[b & 0x0F] ^ hi_tbl[b >> 4]
```

ARM NEON `vqtbl1q_u8` (또는 x86 SSSE3 `pshufb`)가 16바이트를 동시에 테이블 lookup하므로, 한 사이클에 16개 GF(256) 곱셈을 처리할 수 있다.

**gf256_mul_vec 직접 벤치마크** (1MB, 500회 반복, Apple M1 Pro):

| Implementation | Throughput | Speedup |
|---------------|----------:|---------:|
| Scalar (exp/log table) | 1,574 MB/s | 1.0x |
| NEON (vqtbl1q_u8) | 29,674 MB/s | 18.9x |

**인코딩 파이프라인 벤치마크** (RS(5,4), Apple M1 Pro):

| Data size | Scalar (MB/s) | SIMD (MB/s) | Speedup |
|----------|-------------:|------------:|--------:|
| 4 KB | 414 | 4,476 | 10.8x |
| 64 KB | 417 | 6,392 | 15.3x |
| 1 MB | 408 | 6,881 | 16.9x |

소형 데이터에서는 split-table 셋업 오버헤드가 있지만, 64KB 이상에서는 메모리 대역폭에 근접한다. 이 프로젝트에서는 ARM NEON과 x86 SSSE3를 모두 지원하며, 컴파일 시 아키텍처를 자동 감지한다.

참고: 프로덕션 라이브러리(Intel ISA-L 등)는 AVX-512로 64바이트/사이클까지 처리하며, 이 구현은 교육 목적의 단순한 SIMD 적용 예시다.

## Suggested Reading

이 프로젝트는 특정 상용 시스템의 내부 구현을 설명하지 않는다. 배경 지식이 필요하면 공개 문서를 참고하는 편이 안전하다.

- AWS S3 overview: `https://aws.amazon.com/s3/`
- AWS S3 durability docs: `https://docs.aws.amazon.com/AmazonS3/latest/userguide/DataDurability.html`
- AWS S3 performance docs: `https://docs.aws.amazon.com/AmazonS3/latest/userguide/optimizing-performance.html`

## Bottom Line

이 레포는 "S3를 만들었다"는 프로젝트가 아니라, 분산 저장 시스템에서 자주 쓰이는 아이디어 몇 가지를 작은 C 코드로 구현하고 검증한 프로젝트다. 정확성 기준으로 보면 이 방향이 훨씬 정직하고, 엔지니어링 문서로도 더 강하다.
