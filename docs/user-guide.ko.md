# abd-prep 사용자 가이드

Amazon Bio Discovery(ABD, AWS의 항체 설계 서비스)에 올릴 입력 파일을
준비하고 점검하는 명령줄 도구 `abd-prep`의 한국어 상세 가이드입니다.
영어 [README](../README.md)가 기준 문서(레퍼런스)이고, 이 문서는 같은
내용을 **처음 쓰는 사람의 순서대로**, 전문 용어를 풀어 가며 설명합니다.
20분 안에 한 번 돌려 보는 것이 목적이면 더 짧은
[쉬운 사용 방법](easy-tutorial.ko.md)이 있고, 도구를 고치거나 확장하려는
분은 [개발자 가이드](developer-guide.ko.md)를 함께 보세요.

## 목차

- [0. 이 문서를 읽는 방법](#0-이-문서를-읽는-방법)
- [1. abd-prep은 무엇을 하는 도구인가](#1-abd-prep은-무엇을-하는-도구인가)
- [2. 설치](#2-설치)
- [3. 5분 시작: HER2 / 트라스투주맙 예제](#3-5분-시작-her2--트라스투주맙-예제)
- [4. 배경 지식: 파일 안에는 무엇이 들어 있나](#4-배경-지식-파일-안에는-무엇이-들어-있나)
- [5. ABD 레시피별로 필요한 파일](#5-abd-레시피별로-필요한-파일)
- [6. 결과 읽는 법: 메시지 수준, 종료 코드, report.json](#6-결과-읽는-법-메시지-수준-종료-코드-reportjson)
- [7. PDB 구조 파일 다루기 (`abd-prep pdb`)](#7-pdb-구조-파일-다루기-abd-prep-pdb)
- [8. FASTA 서열 파일 다루기 (`abd-prep fasta`)](#8-fasta-서열-파일-다루기-abd-prep-fasta)
- [9. CSV 실험 결과 표 다루기 (`abd-prep csv`)](#9-csv-실험-결과-표-다루기-abd-prep-csv)
- [10. 스크립트와 CI에서 쓰기](#10-스크립트와-ci에서-쓰기)
- [11. 제한 사항](#11-제한-사항)
- [12. 아직 확인되지 않은 ABD 규약 3가지](#12-아직-확인되지-않은-abd-규약-3가지)
- [13. 문제 해결 (FAQ)](#13-문제-해결-faq)
- [14. 용어집](#14-용어집)
- [부록 A. 명령 치트시트](#부록-a-명령-치트시트)
- [부록 B. 자주 만나는 finding code 빠른 대처표](#부록-b-자주-만나는-finding-code-빠른-대처표)

## 0. 이 문서를 읽는 방법

**누구를 위한 문서인가.** ABD에 항체 설계 작업을 올리려는 연구자, 그 일을
도와주는 바이오인포매틱스 담당자, 그리고 실험실에서 나온 결과 표를 ABD
형식으로 정리해야 하는 분입니다. 단백질 구조 파일을 다뤄 본 경험이 없어도
읽을 수 있게 썼습니다. 생물학 용어와 컴퓨터 용어는 처음 나올 때 괄호로
뜻을 적고, 끝의 [용어집](#14-용어집)에 다시 모았습니다.

**명령 표기.** 본문의 명령은 `abd-prep ...` 형태로 적습니다. 설치하지 않고
저장소를 내려받기만 했다면 `abd-prep`을 `python3 -m abd_prep`으로 바꿔
읽으면 됩니다. 두 형태는 완전히 같은 동작을 합니다.

```bash
abd-prep pdb validate my.pdb          # 설치했을 때
python3 -m abd_prep pdb validate my.pdb   # 설치 없이 저장소 디렉터리에서
```

**가장 급한 분을 위한 한 줄 요약.** (터미널이 처음이라면
[쉬운 사용 방법](easy-tutorial.ko.md)부터 보세요.)

| 상황 | 명령 한 줄 |
|---|---|
| ABD를 처음 써 본다 -- 테스트에 필요한 파일을 **한 번에** 전부 만들기 | `abd-prep quickstart` (3.0절) |
| 블로그 예제처럼 1N8Z에서 HER2 타겟 + 트라스투주맙 Fv 파일 만들기 | `abd-prep pdb split --preset her2-trastuzumab --out ./out` |
| 내 PDB 파일이 ABD에 올라갈 상태인지 보기 | `abd-prep pdb validate my.pdb` |
| 내 PDB 파일에서 체인 골라 정리하기 | `abd-prep pdb clean in.pdb out.pdb --keep A:1-107 --keep B:1-120 --rename A=L,B=H` |
| seed 항체 서열(FASTA) 정리하기 | `abd-prep fasta clean seed_raw.fasta seed.fasta --framework igg` |
| 서열 파일과 구조 파일이 같은 항체인지 확인 | `abd-prep fasta validate seed.fasta --framework igg --against fv.pdb` |
| 실험 결과 CSV 점검하기 | `abd-prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit` |
| 실험 결과 CSV를 표준 형태로 저장하기 | `abd-prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit` |

## 1. abd-prep은 무엇을 하는 도구인가

ABD는 세 종류의 파일을 입력으로 받습니다.

| 입력 | 쉽게 말하면 | ABD에서 쓰이는 곳 | 명령 그룹 |
|---|---|---|---|
| **PDB** | 단백질의 3차원 구조. 원자 하나하나의 좌표가 적힌 텍스트 파일 | 설계 목표가 되는 **타겟 항원** 구조, 항체의 **골격(framework) / seed 구조** | `abd-prep pdb` |
| **FASTA** | 단백질의 **서열**. 아미노산을 한 글자씩 이어 쓴 텍스트 파일 | 개량·평가의 출발점이 되는 **seed 항체** | `abd-prep fasta` |
| **CSV** | 엑셀 표를 텍스트로 저장한 것. 한 줄이 항체 변이체 하나, 열이 측정값 | 자기 실험 데이터로 모델을 **미세조정(fine-tuning)** 하거나 결과를 분석할 때 | `abd-prep csv` |

문제는 "뷰어에서는 멀쩡해 보이는데 ABD가 거부하는" 파일이 많다는 점입니다.
한 원자에 좌표가 두 세트 들어 있는 대체 위치(altloc, 결정학자가 원자를 두
자리에서 본 경우), 체인을 지운 뒤 남은 `TER` 레코드, 어느 것이 중쇄이고
어느 것이 경쇄인지 표시가 없는 FASTA, Windows 줄바꿈, 서열 끝의 `*`,
Excel이 몰래 붙이는 UTF-8 BOM(바이트 순서 표식), `;`로 구분된 "CSV" --
이런 것들입니다.

`abd-prep`은 이 **형식 문제**를 고치고, 나머지는 전부 **finding**(발견
사항)이라는 이름의 메시지로 보고합니다. 그리고 한 가지 원칙을 절대로
어기지 않습니다.

> **과학적 내용은 바꾸지 않는다.** 잔기 번호나 이름을 바꾸지 않고, 서열을
> 자르거나 고치지 않고, 측정값을 환산하지 않습니다. 도구가 바꾸는 것은
> 파일의 **포장**(형식·인코딩·구분자·헤더)뿐입니다. 서열에 `X`가 섞여
> 있으면 그것을 지우거나 추측해 바꾸는 대신 "10번째 글자가 X입니다"라고
> 알려 줍니다.

의존성은 없습니다. Python 3.8 이상이면 표준 라이브러리만으로 동작하므로
설치 과정에서 무언가가 깨질 일이 없습니다.

## 2. 설치

### 2.1 설치 없이 바로 실행

```bash
git clone https://github.com/hmkim/abd-prep.git
cd abd-prep
python3 -m abd_prep --help
```

### 2.2 `abd-prep` 명령으로 설치

```bash
pip install .            # 저장소 디렉터리 안에서
abd-prep --version       # abd-prep 0.1.0
```

여러 프로젝트에서 쓰려면 격리 설치 도구 `pipx`가 편합니다.

```bash
pipx install git+https://github.com/hmkim/abd-prep.git
```

### 2.3 네트워크·디스크 사용 범위

- 네트워크를 쓰는 명령은 `pdb fetch`, `pdb split --preset/--pdb`, `quickstart`
  셋뿐입니다.
  RCSB Protein Data Bank(공개 단백질 구조 데이터베이스)의
  `files.rcsb.org`, `data.rcsb.org`에 HTTPS로만 접속합니다.
- 내려받은 파일은 `~/.cache/abd-prep/`에 캐시(한 번 받은 것을 다시 쓰기
  위해 보관)됩니다. 위치는 `--cache DIR` 옵션이나 환경 변수
  `ABD_PREP_CACHE`로 바꿀 수 있습니다.
- 회사 프록시 뒤에 있다면 `HTTPS_PROXY` 환경 변수를 그대로 존중합니다.
- `fasta`·`csv` 그룹은 네트워크를 전혀 쓰지 않습니다.

## 3. 5분 시작: HER2 / 트라스투주맙 예제

ABD 안내 블로그의 §4-3은 PyMOL(구조 뷰어)에서 손으로 하는 절차입니다.
RCSB에서 1N8Z(HER2 세포외 도메인과 허셉틴 Fab의 복합체 결정 구조)를
내려받고, HER2 도메인 IV와 트라스투주맙 Fv만 남기고, 대체 위치와 물을
지우고, 파일 두 개로 저장하는 과정입니다. `abd-prep`은 이것을 한 명령으로
합니다.

### 3.0 가장 간단한 방법: `quickstart`

구조 두 파일뿐 아니라 ABD 테스트에 필요한 **모든** 파일을 한 번에 만들고
싶으면 `quickstart`입니다. 3.1-3.5절의 명령을 순서대로 대신 실행해 주고,
각 파일이 어느 레시피에 쓰이는지를 `README.txt`에 적어 줍니다.

```bash
abd-prep quickstart                      # -> ./abd-quickstart/
abd-prep quickstart --out ./my-test      # 다른 폴더
abd-prep quickstart --no-csv             # CSV 예시 표는 빼고
```

| 만들어지는 파일 | 내용 | 3장에서 해당하는 절 |
|---|---|---|
| `target.pdb` | 타겟 항원 구조(체인 C, 109잔기) | 3.1 (`split`의 `_target.pdb`) |
| `antibody_fv.pdb` | seed 항체 Fv 구조(체인 L·H) | 3.1 (`split`의 `_fv.pdb`) |
| `antibody_seed.fasta` | 같은 항체의 서열(`>H`, `>L`), 구조와 대조까지 끝남 | 3.4 |
| `assay_example.csv` | 결과 표의 **형식 예시**(seed 한 줄 + 자리 표시용 값). 측정값 아님 | 3.5 |
| `README.txt` | 파일별 용도, RCSB 출처·개정일, 다음 단계 | -- |
| `quickstart.report.json` | 전 과정의 기록(출처·SHA-256·발견 사항) | 6.4 |

출력 메시지는 3.1절과 같은 것들에 `seed/`(FASTA 대조)와 `csv/`(표 검증)
접두사가 붙은 줄이 더해진 형태입니다. 같은 폴더에 두 번 실행하면 "이미
파일이 있다"며 종료 코드 2로 멈춥니다(고친 seed 파일을 실수로 덮어쓰지
않도록). 다시 만들려면 `--force`. 프리셋은 `--preset`으로 고르고, 자기
구조·서열·표에는 아래 3.6절과 7-9장의 개별 명령을 씁니다. 실제 출력
전문과 한 줄씩 설명은 [쉬운 사용 방법 4장](easy-tutorial.ko.md#4-한-번에-준비하기-quickstart-1분)에
있습니다.

### 3.1 한 명령으로 잘라내고, 정리하고, 검증하기

```bash
abd-prep pdb split --preset her2-trastuzumab --out ./out
```

실제 출력입니다.

```
[split] her2-trastuzumab  ->  target: ./out/her2-trastuzumab_target.pdb, fv: ./out/her2-trastuzumab_fv.pdb
  WARNING PDB_RESOLUTION       1N8Z              resolution 2.52 A is coarser than the 2.5 A the ABD guide recommends
  INFO    PDB_HETATM_REMOVED   target/file       removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    PDB_GAP_NUMBERING    target/C:580-591  residue numbering jumps 580 -> 591 (expected)
  INFO    PDB_GAP_CA_DISTANCE  target/C:580-591  CA-CA distance 6.7 A between consecutive residues (> 4.2 A): chain break -- expected for this preset
  INFO    PDB_HETATM_REMOVED   fv/file           removed HETATM records: HOH x79, NAG x28, SO4 x5
  summary: 0 error(s), 1 warning(s), 4 info
  report: ./out/her2-trastuzumab.report.json
```

한 줄씩 읽어 보겠습니다.

- 첫 줄: 프리셋(미리 정해 둔 선택 규칙) `her2-trastuzumab`으로 파일 두 개를
  만들었습니다. `_target.pdb`가 **항원**, `_fv.pdb`가 **항체 Fv**입니다.
- `WARNING PDB_RESOLUTION`: 1N8Z의 해상도(구조가 얼마나 선명하게 측정됐는지,
  숫자가 작을수록 좋음)는 2.52 Å인데 ABD 가이드는 2.5 Å 이하를 권장합니다.
  사실이므로 그대로 보고합니다. 경고(WARNING)는 파일을 못 쓰게 하지는
  않습니다.
- `INFO PDB_HETATM_REMOVED`: 물(HOH) 79개, 글리칸 NAG(단백질에 붙은 당)
  28개, 황산 이온(SO4) 5개를 지웠습니다. 정보(INFO)는 "이런 일을 했다"는
  기록입니다.
- `INFO PDB_GAP_NUMBERING` / `PDB_GAP_CA_DISTANCE`: 타겟 체인 C의 잔기
  580번 다음이 591번입니다. 1N8Z 결정에서 581-590번 고리(loop)가 보이지
  않아 모델에 빠져 있기 때문이고, 프리셋이 "이 자리는 원래 비어 있다"고
  알고 있어서 `(expected)` 정보로만 표시합니다. 프리셋 없이 돌리면 같은
  자리가 경고로 나옵니다(7.5절).
- `summary`: 오류 0, 경고 1, 정보 4. 종료 코드 0(통과).
- `report`: 무엇을 했는지, 어느 RCSB 개정판을 썼는지, 각 파일의 SHA-256
  (파일 내용의 지문)이 JSON으로 남습니다(6.4절).

### 3.2 만들어진 파일 들여다보기

```bash
abd-prep pdb inspect ./out/her2-trastuzumab_fv.pdb
```

```
== ./out/her2-trastuzumab_fv.pdb  altloc atoms: 0  HETATM: none  MODEL: 0  TER: 2  END: yes
  chain L: 1-107 (107 res, 823 atoms)
     DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
  chain H: 1-120 (120 res, 924 atoms)
     EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
```

대체 위치 0, HETATM 없음, 체인 경계마다 `TER` 하나(2개), `END` 있음.
경쇄 L 107잔기, 중쇄 H 120잔기 -- 트라스투주맙 Fv의 문헌값과 일치합니다.
이 두 서열은 구조에서 바로 읽어 낸 것이라, 다음 단계의 FASTA 출발점으로
쓸 수 있습니다.

### 3.3 업로드

- `her2-trastuzumab_target.pdb` -> ABD의 **타겟 항원** 파일(항원을 받는
  모든 레시피).
- `her2-trastuzumab_fv.pdb` -> **항체 골격 / seed 구조** 파일(유도 진화와
  속성 예측은 seed를 PDB 또는 FASTA로 받습니다).

### 3.4 seed를 서열(FASTA)로도 준비하기

seed를 서열로 올리고 싶다면 구조에서 뽑아낸 뒤, 두 파일이 같은 분자인지
확인까지 해 둡니다.

```bash
abd-prep pdb inspect ./out/her2-trastuzumab_fv.pdb --fasta > seed.fasta
abd-prep fasta validate seed.fasta --framework igg --against ./out/her2-trastuzumab_fv.pdb
```

```
[fasta validate] seed.fasta
  INFO    FASTA_PDB_IDENTITY  H:H<-H  PDB chain H vs record 'H': 120/120 aligned residues identical (PDB residue 1 = FASTA position 1)
  INFO    FASTA_PDB_IDENTITY  L:L<-L  PDB chain L vs record 'L': 107/107 aligned residues identical (PDB residue 1 = FASTA position 1)
  summary: 0 error(s), 0 warning(s), 2 info
```

논문이나 데이터베이스에서 받은 서열이 따로 있다면 그 파일을 정리합니다
(8장).

```bash
abd-prep fasta clean trastuzumab_from_paper.fasta seed.fasta --framework igg
```

### 3.5 내 실험 결과 표

```bash
abd-prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv validate results_clean.csv --against seed.fasta    # 표의 행들이 seed를 가리키는가?
```

### 3.6 프리셋 없이 직접 지정하기

다른 RCSB 항목이나 내 파일이라면 체인과 잔기 범위를 직접 씁니다.

```bash
abd-prep pdb split --pdb 1N8Z --target C:489-607 --heavy B:1-120 --light A:1-107 \
    --rename-heavy H --rename-light L --out ./out
abd-prep pdb split --input ./my_complex.pdb --target B:530-601 --out ./out
```

## 4. 배경 지식: 파일 안에는 무엇이 들어 있나

이 절은 세 파일 형식을 처음 보는 분을 위한 10분짜리 설명입니다. 이미
익숙하면 5장으로 넘어가세요.

### 4.1 PDB 파일 -- 원자 좌표 목록

PDB 파일은 1970년대부터 쓰인 고정폭 텍스트 형식입니다. 한 줄이 원자
하나입니다.

```
ATOM   1667  N   ASP L   1      37.123  12.456  -3.789  1.00 20.50           N
ATOM   1668  CA  ASP L   1      38.001  11.900  -4.800  1.00 21.10           C
HETATM 9001  O   HOH A 301      10.000  10.000  10.000  1.00 30.00           O
TER    9002      LYS L 107
END
```

- **`ATOM`**: 단백질(또는 핵산)의 원자. **`HETATM`**: 그 외 모든 것 --
  물(HOH), 이온(SO4), 당(NAG), 약물 분자, 그리고 **변형된 아미노산**
  (MSE=셀레노메티오닌, CSO=산화 시스테인 등). 마지막 종류가 까다로운데,
  단백질 사슬의 일부인데도 `HETATM`으로 적혀 있어 "HETATM을 다 지우라"는
  규칙을 그대로 따르면 사슬에 구멍이 납니다(7.3절).
- **체인(chain)**: 한 폴리펩타이드 사슬. 22번째 칸의 글자 하나(`L`, `H`,
  `A`...)로 구분합니다. 항체는 중쇄(heavy)와 경쇄(light) 두 체인, 항원은
  보통 한 체인입니다.
- **잔기(residue)**: 사슬을 이루는 아미노산 하나. 세 글자 이름(`ASP`)과
  **잔기 번호**(23-26번째 칸)가 있고, 번호 뒤에 **삽입 코드**(27번째 칸,
  `82A`의 `A`)가 붙을 수 있습니다. 항체 분야에서 많이 쓰는 Kabat 번호
  체계는 이 삽입 코드로 `82`, `82A`, `82B`, `82C`처럼 번호를 매깁니다.
- **대체 위치(altloc, alternate location)**: 17번째 칸의 글자. 결정학자가
  한 원자를 두 자리에서 봤을 때 `A`/`B`로 두 줄을 넣습니다. 각 줄의
  **점유율(occupancy, 55-60번째 칸)** 은 그 자리에 있을 확률(0.6/0.4처럼
  합이 1)입니다. ABD는 이 상태의 파일을 "PDB contains alternate
  locations"라며 거부합니다. 정리는 둘 중 하나를 고르는 일이고, `abd-prep`은
  **점유율이 높은 쪽**을 고릅니다(PyMOL의 "A를 남기기"와 다를 수 있음).
- **`TER`**: 체인이 끝났다는 표시. 체인을 지웠는데 그 체인의 `TER`이 남아
  있으면 엄격한 파서가 걸려 넘어집니다. **`END`**: 파일 끝.
- **`MODEL`/`ENDMDL`**: NMR 구조처럼 모델이 여러 개 들어 있을 때. ABD
  입력은 모델 하나여야 하므로 첫 모델만 남깁니다.
- **해상도(resolution, Å=옹스트롬, 1Å = 0.1nm)**: 구조가 얼마나 선명한지.
  작을수록 좋고 ABD 가이드는 2.5 Å 이하를 권장합니다.
- **Fab / Fv / VH / VL**: 항체(Y자 모양)에서 항원을 붙잡는 끝부분이
  Fab, 그 중에서도 실제 결합을 담당하는 가변 영역 두 도메인(중쇄 쪽
  **VH**, 경쇄 쪽 **VL**)을 합쳐 **Fv**라고 합니다. ABD의 항체 파일은 이
  Fv만 담는 것이 보통이고, 공개 결정 구조에는 Fab 전체(불변영역 포함)가
  들어 있어 잘라내야 합니다.

### 4.2 FASTA 파일 -- 서열을 글자로

```
>H heavy trastuzumab VH
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG...
>L VL_1 kappa trastuzumab VL
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSG...
```

- `>`로 시작하는 줄이 **헤더**, 그 아래 줄들이 **서열**입니다. 헤더 +
  서열 한 묶음이 **레코드**(record)이고, 헤더의 첫 단어가 레코드
  **ID**입니다.
- 서열은 아미노산 20종의 **1문자 기호**(A, C, D, E, F, G, H, I, K, L, M,
  N, P, Q, R, S, T, V, W, Y)로만 써야 합니다. `X`(모름), `B`, `Z`, `J`,
  `U`, `O`는 모호한 글자로, ABD가 "non-standard residue"로 거부합니다.
- 끝의 `*`는 DNA를 번역할 때 생기는 **종결 코돈** 표시입니다. 끝에 있으면
  지워도 되지만, **중간에** 있으면 번역이 잘못된 서열이라는 뜻이므로
  지우면 안 됩니다.
- 항체 FASTA에서 가장 흔한 사고는 **중쇄와 경쇄를 바꿔 올리는 것**입니다.
  그래서 `abd-prep`은 헤더에 분명한 표시가 없으면 절대 추측하지 않고
  사용자에게 지정하라고 요구합니다(8.2절).

### 4.3 CSV 파일 -- 실험 결과 표

```
id,sequence,value,unit,target,assay,notes
Ab001,EVQLVESGGG...,2.3,nM,HER2,SPR,baseline
Ab002,EVQLVESGGG...,0.8,nM,HER2,SPR,
```

- 쉼표(또는 `;`, 탭)로 열을 나눈 텍스트. 첫 줄이 **헤더**(열 이름).
- **인코딩**: 글자를 바이트로 적는 방식. ABD는 **UTF-8**만 받습니다.
  한국어 Windows의 Excel이 "CSV"로 저장하면 CP949(EUC-KR)로 저장되기 쉽고,
  "CSV UTF-8"로 저장하면 앞에 **BOM**(바이트 순서 표식, 보이지 않는 3바이트)
  이 붙습니다. 둘 다 `abd-prep csv clean`이 처리합니다.
- **구분자**: 한국어·유럽 로캘의 Excel은 `;`로 구분해 저장하는 경우가
  많습니다. 역시 `clean`이 `,`로 바꿉니다.
- 측정값 예: **KD**(결합 세기, 작을수록 강함, nM·pM 단위), **Tm**(열
  안정성, degC), 발현량(mg/L). 장비 소프트웨어(SPR, BLI, ELISA 판독기)가
  내보내거나 수탁기관(CRO)이 Excel로 보내 줍니다.

## 5. ABD 레시피별로 필요한 파일

ABD는 작업 흐름을 **레시피**(recipe)로 묶습니다. 어떤 파일이 필요한지는
"이미 있는 항체(**seed**)에서 출발하는가"로 갈립니다.

| 레시피 | 하는 일 | 타겟 항원 | seed 항체 | 비고 |
|---|---|---|---|---|
| **Antibody generation -- De novo** | 처음부터 새 항체 생성 | **필수**(PDB) | 안 씀 | 골격(framework) 선택: Nanobody / IgG / scFv |
| **Antibody generation -- Directed evolution** | 기존 항체 개량 | 선택(PDB) | **필수**(FASTA 또는 PDB) | 항원을 넣으면 결합력 중심 최적화, 빼면 항체 자체 속성 최적화 |
| **Property prediction** | 기존 항체의 결합력·안정성·개발 적합성 예측 | 선택 | **필수**(FASTA 또는 PDB) | |

**골격(framework)** 은 seed 파일에 체인이 몇 개 있어야 하는지를
결정하고, `abd-prep fasta`는 `--framework` 옵션으로 이를 검사합니다.

| framework | 무엇인가 | seed 파일의 레코드 |
|---|---|---|
| **IgG** | 사람 표준 Y자 항체. 가변 영역이 중쇄 VH + 경쇄 VL 두 도메인 | **heavy 1개 + light 1개** |
| **Nanobody (VHH)** | 낙타류의 중쇄-only 항체에서 온 **단일 도메인** 항체 | **1개**(heavy만) |
| **scFv** | VH와 VL을 유연한 **링커**(보통 `GGGGS` 반복)로 이어 붙인 인공 단일 사슬 | **1개**(두 도메인과 링커가 한 서열에) |

모든 업로드 파일은 **100 MB 미만**, **UTF-8**이어야 하며, `abd-prep`은 어떤
명령이든 이 둘을 먼저 확인합니다.

## 6. 결과 읽는 법: 메시지 수준, 종료 코드, report.json

### 6.1 메시지 수준 세 가지

모든 명령은 같은 모양의 표를 출력합니다. 각 줄은 `수준 CODE 위치 설명`
입니다.

| 수준 | 뜻 | 종료 코드에 미치는 영향 |
|---|---|---|
| `ERROR` | ABD가 거부할 것이 거의 확실한 문제, 또는 도구가 고칠 수 없는 문제 | 1 |
| `WARNING` | 확인해 볼 만한 점. 파일이 틀렸다는 뜻은 아님 | 0 (`--strict`를 주면 1) |
| `INFO` | 도구가 한 일, 또는 참고 사항 | 없음 |

`CODE`(예: `PDB_ALTLOC_PRESENT`)는 **finding code**라고 부르는 고정
식별자입니다. 메시지 문구는 바뀔 수 있지만 코드는 바뀌지 않으므로,
스크립트에서 걸러 내거나 README의 [Finding codes](../README.md#finding-codes)
표에서 찾을 때 이 코드를 쓰세요. 코드는 `FILE_`(모든 형식 공통), `PDB_`,
`FASTA_`, `CSV_` 접두로 시작합니다.

### 6.2 종료 코드

| 코드 | 뜻 |
|---|---|
| **0** | 통과. 경고와 정보는 있어도 됨 |
| **1** | 검증 오류가 하나 이상 (또는 `--strict`에서 경고가 하나 이상) |
| **2** | 입력·옵션 오류: 파일 없음, 모르는 체인, 잘못된 범위, CSV 자리에 `.xlsx` 등. 메시지는 `abd-prep: error: ...` 한 줄로 표준 오류에 출력 |

`clean`은 검증 오류가 있어도 **출력 파일을 씁니다**. 결과를 눈으로 보고
판단할 수 있게 하려는 것이고, 종료 코드는 그대로 1입니다. 유일한 예외는
파일을 아예 읽어 낼 수 없는 경우(인코딩 오류)입니다.

### 6.3 모든 명령이 받는 공통 옵션

| 옵션 | 뜻 |
|---|---|
| `--json` | 텍스트 표 대신 리포트 전체를 JSON으로 출력 |
| `--strict` | 경고도 실패(종료 1)로 취급. CI·스크립트용 |
| `--encoding ENC` | FASTA/CSV 입력의 인코딩(기본 `utf-8`; `cp949`, `latin-1` 등). `clean`은 UTF-8로 변환해 저장 |
| `--report PATH` | `report.json`을 PATH에 저장 |
| `--no-report` | `report.json`을 쓰지 않음 (`clean`은 기본으로 `<출력 파일 이름>.report.json`을 씁니다. `validate`는 `--report`를 줄 때만 씀) |
| `--report-paths basename\|full` | 리포트에 파일 경로를 파일 이름만(기본) 적을지 절대 경로로 적을지. 리포트는 메일·티켓에 첨부되기 쉬워 기본값이 파일 이름입니다 |

### 6.4 report.json -- 재현성 기록

`clean`과 `split`은 출력 옆에 리포트를 남기고, `validate`는 `--report`를
주면 남깁니다. 세 형식 모두 같은 틀입니다. 3.4절의 `fasta clean`이 남긴
실제 리포트(앞부분)입니다.

```json
{
  "tool": {"name": "abd-prep", "version": "0.1.0", "python": "3.9.25"},
  "command": "abd-prep fasta clean seed_igg.fasta seed.fasta --framework igg",
  "generated_at": "2026-10-26T19:19:35Z",
  "inputs":  [{"kind": "fasta", "path": "seed_igg.fasta", "sha256": "85ee1dec...", "bytes": 478}],
  "outputs": [{"kind": "fasta", "path": "seed.fasta",     "sha256": "7970e6fe...", "bytes": 471}],
  "findings": [
    {"level": "info", "code": "FASTA_GAPS_REMOVED",  "message": "alignment gap character(s) '-'/'.' removed (1 run(s))", "where": "heavy"},
    {"level": "info", "code": "FASTA_LOWERCASE",     "message": "14 lower-case letter(s) upper-cased", "where": "heavy"},
    {"level": "info", "code": "FASTA_TRAILING_STOP", "message": "terminal '*' (translation stop) removed", "where": "VL_1"},
    {"level": "info", "code": "FASTA_CRLF",          "message": "6 line(s) end with CR/LF (Windows); written as LF", "where": "file"}
  ],
  "details": {"framework": "igg", "records": [{"id": "heavy", "length": 120, "role": "H", "role_source": "header", "lineno": 1}, "..."]}
}
```

- `command`: 실제로 실행한 명령. 몇 달 뒤 "이 파일을 어떻게 만들었지?"에
  답합니다.
- `inputs`/`outputs`의 `sha256`: 파일 내용의 지문. 바이트 하나라도 다르면
  값이 달라지므로, ABD에 올린 파일이 이 리포트의 파일과 같은지 확인할 수
  있습니다.
- `details`: 명령별 세부 정보. `split`은 RCSB 항목 정보(개정일, 해상도,
  실험 방법)와 파트별(target/fv) 정리·검증 결과를, `fasta`는 레코드 목록과
  `--against` 비교 결과를, `csv`는 열 매핑·결측 개수·단위·행 수를 담습니다.

이 리포트를 결과 파일과 함께 보관해 두면 실험 기록(ELN)에 "입력 출처"
항목으로 그대로 붙일 수 있습니다.

## 7. PDB 구조 파일 다루기 (`abd-prep pdb`)

```
abd-prep pdb fetch    1N8Z 6LBX [--out DIR] [--cache DIR] [--refresh] [--no-meta]
abd-prep pdb split    --preset NAME | --pdb ID | --input FILE  [--target SPEC] [--heavy SPEC] [--light SPEC]
                      [--rename-target ID] [--rename-heavy ID] [--rename-light ID] [--name NAME] [--out DIR]
                      [--keep-het] [--std-residues] [--strip-h] [--renumber]
abd-prep pdb clean    IN.pdb OUT.pdb [--keep SPEC ...] [--rename A=L,B=H] [--keep-het] [--std-residues] [--strip-h] [--renumber]
abd-prep pdb validate FILE... [--preset NAME] [--expect-chains N] [--allow-het]
abd-prep pdb inspect  FILE... [--fasta] [--json]
```

### 7.1 상황별로 고르는 명령

| 하고 싶은 일 | 명령 |
|---|---|
| 공개 구조(RCSB)에서 타겟 + 항체 Fv 두 파일을 한 번에 만들기 | `pdb split --pdb ID --target ... --heavy ... --light ...` (프리셋이 있으면 `--preset NAME`) |
| 내 PDB 파일(모델링 결과, 자체 결정 구조)에서 두 파일 만들기 | `pdb split --input FILE --target ... --heavy ... --light ...` |
| 파일 하나를 세밀하게 정리하기 (체인 고르기, 이름 바꾸기, 옵션 조합) | `pdb clean IN OUT --keep ... --rename ...` |
| 업로드 전에 점검만 하기 (파일을 바꾸지 않음) | `pdb validate FILE` |
| 체인 구성·잔기 범위·서열을 들여다보기 | `pdb inspect FILE` (`--fasta`로 서열을 FASTA로 뽑기) |
| RCSB 파일을 내려받아 두기 | `pdb fetch ID ...` |

`split`은 `fetch` + `clean` + `validate`를 합친 것입니다. 결과는
`<이름>_target.pdb`(체인 1개), `<이름>_fv.pdb`(heavy + light),
`<이름>.report.json`이고, 이름은 프리셋 이름, 항목 ID, 또는 입력 파일 이름에서
자동으로 정해집니다(`--name`으로 지정 가능). 타겟과 항체를 서로 다른 항목에서
가져올 수도 있습니다(프리셋 `her2-trastuzumab-6lbx-6b9z`가 그 예).

### 7.2 체인과 잔기 고르는 문법

`--keep`, `--target`, `--heavy`, `--light`는 모두 같은 문법을 씁니다.

| 쓰는 법 | 뜻 |
|---|---|
| `B` | 체인 B 전체 |
| `A:1-107` | 체인 A의 잔기 1번부터 107번까지 |
| `C:489-607,610-612` | 범위 여러 개 |
| `C:605` | 잔기 하나 |
| `--keep A:1-107 --keep B:1-120` | (`clean`) 체인 여러 개는 옵션을 반복 |

번호는 **파일에 적힌 그대로의 잔기 번호**입니다(순차 번호가 아님). 삽입
코드가 있는 Kabat 번호 체계도 그대로 다룹니다.

- 상한을 삽입 코드 **없이** 쓰면 그 번호의 삽입 잔기를 모두 포함합니다:
  `H:1-82`는 82, 82A, 82B, 82C를 모두 남깁니다.
- 상한에 삽입 코드를 **붙여** 쓰면 정확히 거기까지입니다: `H:1-82A`는
  82와 82A까지.
- 음수 번호(`-3`)도 받습니다.

### 7.3 `clean`이 하는 일, 순서대로

1. **첫 `MODEL`만 남깁니다.** 모델이 여럿이면 경고 `PDB_MULTI_MODEL`.
2. **체인·잔기 선택**을 적용합니다. 없는 체인을 지정하면
   `PDB_CHAIN_NOT_FOUND`, 범위가 아무것도 고르지 않으면
   `PDB_SELECTION_EMPTY` 오류.
3. **대체 위치를 정리합니다.** 원자마다 점유율이 가장 높은 컨포머(형태)를
   남기고 17번째 칸을 비웁니다. 동률이면 파일에서 먼저 나온 것. 한 자리에
   **서로 다른 잔기 종류**가 모델링된 경우(예: 6B9Z B175의 `ACSO`/`BCYS`)는
   점유율이 높은 잔기 이름을 먼저 고르고, HETATM 필터에서 지워질 컨포머가
   살아남을 컨포머를 밀어내지 않게 합니다(`PDB_RESIDUE_HETEROGENEITY` 정보).
4. **HETATM을 지웁니다**(`--keep-het`를 주면 남김). 물·이온·당·리간드가
   여기서 사라집니다. 단, 단백질 백본의 일부인 변형 잔기(N, CA, C 원자를
   갖고 사슬 안에 있는 것 -- MSE, CSO, PCA ...)를 지우면 사슬에 구멍이 나므로
   **경고** `PDB_MODIFIED_RESIDUE_REMOVED`를 냅니다. 구멍을 내지 않으려면
   `--std-residues`로 표준 아미노산으로 **변환**합니다.

   | 변형 잔기 | 변환 결과 | 원자 처리 |
   |---|---|---|
   | `MSE`(셀레노메티오닌) | `MET` | `SE` -> `SD`로 이름 변경, 원소 칸 갱신 |
   | `CSO`, `CSD`, `OCS`, `CME`(산화·수식 시스테인) | `CYS` | 산소·부가 원자 삭제 |
   | `SEP`, `TPO`, `PTR`(인산화) | `SER`, `THR`, `TYR` | 인산기 삭제 |
   | `HYP`(하이드록시프롤린) | `PRO` | `OD1` 삭제 |
   | `MLY`, `KCX`(수식 라이신) | `LYS` | 메틸·카복실 원자 삭제 |
   | `PCA`(피로글루탐산) | `GLN` | `OE` -> `OE1`; 측쇄는 불완전한 채로 남음 |

   표준 잔기의 원자 집합에 없는 원자는 지우고 리포트에 적습니다. 변환
   규칙이 없는 변형 잔기는 경고(`PDB_MODIFIED_RESIDUE_NO_RULE`)와 함께
   지워지며 `--keep-het`로 남길 수 있습니다.
5. **수소 제거**(`--strip-h`): 수소·중수소 원자를 원소 칸 기준으로(칸이
   비어 있으면 원자 이름으로) 지웁니다.
6. **체인 이름 변경**(`--rename A=L,B=H`): 쓰는 시점에 줄 단위로 바꾸므로
   두 ID를 서로 맞바꾸는 것(`A=B,B=A`)도 안전합니다. 두 체인을 한 ID로
   합치려 하면 `PDB_RENAME_COLLISION` 오류.
7. **`TER`과 `END`를 다시 만듭니다.** 체인 경계마다 마지막 잔기를 적은
   `TER` 하나. 원본의 `TER`은 전부 버립니다(지운 체인의 `TER`이 남는 것이
   엄격한 파서가 거부하는 흔한 이유). 원자 일련번호는 원본 그대로 두고
   (`--renumber`를 주면 1부터 다시 매김), 중복이 생기면 `validate`가
   `PDB_SERIAL_DUPLICATE` 경고.

이 외의 `ATOM`/`HETATM` 내용(좌표, B-factor, 원자 순서)은 **바이트 그대로**
복사됩니다. `HEADER`, `REMARK`, `SEQRES`, `SSBOND`, `LINK`, `CRYST1` 등
다른 레코드는 쓰지 않습니다(11장).

**실제 예 -- 원본 1N8Z에서 Fv 잘라내기**

```bash
abd-prep pdb fetch 1N8Z --out .
abd-prep pdb clean 1N8Z.pdb her2_fv.pdb --keep A:1-107 --keep B:1-120 --rename A=L,B=H
```

```
[clean] 1N8Z.pdb -> her2_fv.pdb  chains: L(107 res), H(120 res)
  INFO    PDB_HETATM_REMOVED  file  removed HETATM records: HOH x79, NAG x28, SO4 x5
  summary: 0 error(s), 0 warning(s), 1 info
  report: her2_fv.report.json
```

```bash
abd-prep pdb validate her2_fv.pdb --expect-chains 2
```

```
[validate] her2_fv.pdb
  (no findings)
  summary: 0 error(s), 0 warning(s), 0 info
```

### 7.4 프리셋

프리셋은 "어느 항목의 어느 체인·범위를 어떤 이름으로 쓸지"와 "결과가
어떠해야 하는지"(잔기 수, 알려진 결손)를 묶어 둔 것입니다.

| 프리셋 | 항목 | 타겟 | Fv | 비고 |
|---|---|---|---|---|
| `her2-trastuzumab` (기본) | 1N8Z (HER2 ECD / 허셉틴 Fab, 2.52 Å) | `C:489-607` -> 체인 `C`, 109잔기 | `A:1-107` -> `L`, `B:1-120` -> `H` | 블로그 절차·트라스투주맙 표준 서열과 일치. 체인 C는 성숙 단백질 번호(UniProt P04626 번호 - 22). 581-590 고리는 결정에 없어 "예상된 결손"으로 처리 |
| `her2-trastuzumab-6lbx-6b9z` | 6LBX(타겟) + 6B9Z(Fv) | `B:530-601` -> 체인 `B`, 72잔기 | `A:1-107` -> `L`, `B:1-120` -> `H` | 6B9Z 항체는 트라스투주맙과 프레임워크 몇 자리가 달라 "트라스투주맙 유래 변이체"로 볼 것. 타겟과 항체가 다른 결정에서 왔으므로 좌표계는 서로 무관 |

- `validate --preset NAME`은 파일의 체인 ID 집합을 프리셋의 target/fv 파트와
  맞춰 보고 기대치(잔기 수·예상 결손)를 적용합니다. `split` 밖에서 만든
  파일도 같은 기준으로 검사할 수 있습니다.

  ```bash
  abd-prep pdb validate ./out/her2-trastuzumab_target.pdb ./out/her2-trastuzumab_fv.pdb --preset her2-trastuzumab
  ```

  ```
  [validate] ./out/her2-trastuzumab_target.pdb ./out/her2-trastuzumab_fv.pdb
    INFO    PDB_GAP_NUMBERING    her2-trastuzumab_target.pdb:C:580-591  residue numbering jumps 580 -> 591 (expected)
    INFO    PDB_GAP_CA_DISTANCE  her2-trastuzumab_target.pdb:C:580-591  CA-CA distance 6.7 A between consecutive residues (> 4.2 A): chain break -- expected for this preset
    summary: 0 error(s), 0 warning(s), 2 info
  ```

- 두 프리셋 모두 타겟 체인은 **원본 ID를 유지**합니다(ABD가 타겟 체인을
  `T`로 요구하는지 아직 확인되지 않음, 12장). 바꾸려면 `--rename-target T`.
- **내 항체의 Fv 경계는?** 프리셋의 107/120은 트라스투주맙의 문헌값입니다.
  다른 항체는 `--heavy`/`--light`에 범위를 직접 넣어야 하고, 도구는 Fv
  경계를 계산해 주지 않습니다(ANARCI 같은 항체 번호 매기기 도구는 포함되지
  않음). 모르면 `pdb inspect`로 서열을 보고 VH는 보통 `...WGQGT(L/T)VTVSS`,
  VL은 `...FGQGTK(V/L)EIK` 근처에서 끝난다는 점을 참고하세요.

### 7.5 `validate`가 보는 것과 자주 보는 메시지

`validate`는 ABD 파서가 걸려 넘어지는 것들을 봅니다: 대체 위치, HETATM,
`TER`/`END` 위치, 체인 수, VL(105-115)·VH(110-125)로서 그럴듯한 잔기 수,
사슬 끊어짐(연속 잔기의 CA-CA 거리 > 4.2 Å; 정상 펩타이드 결합은 약
3.8 Å), 잔기 번호 점프, 비표준 잔기 이름, 중복 일련번호, 비ASCII 바이트.

**원본 RCSB 파일을 그대로 검사하면** 이렇게 나옵니다(1N8Z).

```
[validate] 1N8Z.pdb
  WARNING PDB_HETATM_PRESENT   file       HETATM records present: HOH x79, NAG x28, SO4 x5; ABD inputs normally contain protein atoms only (clean removes them by default)
  INFO    PDB_GAP_NUMBERING    C:101-111  residue numbering jumps 101 -> 111
  WARNING PDB_GAP_CA_DISTANCE  C:101-111  CA-CA distance 15.2 A between consecutive residues (> 4.2 A): chain break
  ...
  ERROR   PDB_TER_MISSING      A          no TER after residue 233
  ERROR   PDB_TER_MISPLACED    line 8486  TER (ASN C 607) is not at a chain boundary -- left over from a removed chain or cut residue
  summary: 5 error(s), 5 warning(s), 4 info
```

원본 파일은 각 체인의 `TER` 뒤에 그 체인 ID를 단 물 분자(HETATM)가
이어지므로 `TER` 규칙에 걸립니다. 이것은 "원본이 틀렸다"가 아니라
"업로드용 파일이 아니다"라는 뜻이고, `clean`을 거치면 전부 사라집니다
(7.3절의 `(no findings)`). `validate`는 **업로드 직전 파일**에 쓰는
명령이라고 생각하면 됩니다.

| 메시지 | 뜻 | 대처 |
|---|---|---|
| `ERROR PDB_ALTLOC_PRESENT` | 대체 위치가 남아 있음. ABD의 "contains alternate locations" 거부 그 자체 | `pdb clean` |
| `WARNING PDB_HETATM_PRESENT` | 물·이온·당 등이 남아 있음 | `pdb clean` (의도적으로 남긴 것이면 `validate --allow-het`) |
| `ERROR PDB_TER_MISSING` / `PDB_TER_MISPLACED` / `PDB_TER_EXTRA` / `PDB_TER_MISMATCH` | `TER`이 없거나, 지운 체인의 것이 남았거나, 겹치거나, 앞 잔기와 다른 잔기를 가리킴 | `pdb clean`이 `TER`을 다시 만듦 |
| `WARNING PDB_GAP_CA_DISTANCE` | 연속 잔기 사이가 멀다 = 결정에서 보이지 않은 고리가 빠져 있음 | 사실이면 그대로 둠. ABD는 결손이 있는 구조도 받지만 그 자리의 설계는 신뢰가 낮음. 프리셋이 알고 있는 결손은 정보로만 표시 |
| `INFO PDB_GAP_NUMBERING` | 잔기 번호가 뛰거나 줄어듦(Kabat 삽입 코드, 결손, 키메라) | 참고용 |
| `WARNING PDB_RESIDUE_COUNT` | 체인 H/L의 잔기 수가 Fv로서 어색함 | 범위가 너무 짧거나(잘못 잘림) 길지(불변영역 포함) 않은지 확인 |
| `WARNING PDB_MODIFIED_RESIDUE_REMOVED` | 선택 범위 안의 변형 잔기가 HETATM과 함께 지워져 구멍이 생김 | `clean --std-residues` |
| `WARNING PDB_NONSTANDARD_RESNAME` | `ATOM` 레코드에 표준 아미노산이 아닌 이름 | 모델링 도구의 산출물인지 확인 |
| `WARNING PDB_RESOLUTION` (`split`) | 항목 해상도가 2.5 Å보다 거침 | 사실 보고. 더 좋은 구조가 있는지 검토 |
| `WARNING PDB_SERIAL_DUPLICATE` | 원자/`TER` 일련번호 중복 | `clean --renumber` |
| `WARNING PDB_MULTI_MODEL` | 모델이 여럿 | `clean`이 첫 모델만 남김 |
| `WARNING PDB_NON_ASCII` | `ATOM` 줄에 비ASCII 바이트 | 편집기에서 유니코드 문자가 섞인 경우. 줄을 확인 |

### 7.6 `inspect`와 `fetch`

`pdb inspect`는 파일을 바꾸지 않고 요약합니다: 대체 위치 원자 수, HETATM
종류, 모델 수, `TER`/`END`, 체인별 잔기 범위·개수·원자 수·서열.
`--fasta`를 주면 체인 서열을 FASTA로 출력하는데, 헤더에 ABD 규약의 `>H`/`>L`
(체인 ID가 H/L일 때)을 쓰므로 그대로 `fasta validate`에 넣을 수 있습니다.

```bash
abd-prep pdb inspect ./out/her2-trastuzumab_fv.pdb --fasta
```

```
>L her2-trastuzumab_fv.pdb chain L 1-107 (107 res)
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
>H her2-trastuzumab_fv.pdb chain H 1-120 (120 res)
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
```

`pdb fetch`는 RCSB에서 PDB 형식 파일을 받아 캐시에 두고(이미 있으면 다시
받지 않음), `--out DIR`로 사본을 복사하고, 항목 정보(제목, 실험 방법,
해상도, 공개일, 개정일)를 함께 보여 줍니다.

```
1N8Z: ~/.cache/abd-prep/1N8Z.pdb  (703809 bytes, cache)
  Crystal structure of extracellular domain of human HER2 complexed with Herceptin Fab
  method X-ray, resolution 2.52 A, released 2003-02-18, revised 2024-10-30
```

RCSB가 mmCIF 형식으로만 제공하는 아주 큰 항목(체인 62개 또는 원자 99,999개
초과)은 "mmCIF 전용"이라는 메시지로, 존재하지 않는 ID는 "없는 항목"으로
구분해 알려 줍니다. 항목 정보를 못 받아도(네트워크 제한 등) 파일 처리는
그대로 진행하고 `PDB_ENTRY_META_UNAVAILABLE` 정보만 남깁니다.

## 8. FASTA 서열 파일 다루기 (`abd-prep fasta`)

```
abd-prep fasta clean    IN.fasta OUT.fasta --framework igg|nanobody|scfv [--heavy-id ID] [--light-id ID]
                        [--replace-ambiguous X=G,...] [--no-heuristics] [--wrap N]
abd-prep fasta validate FILE... [--framework ...] [--heavy-id ID] [--light-id ID] [--replace-ambiguous ...]
                        [--no-heuristics] [--against FV.pdb] [--chain-map H=B,L=A]
abd-prep fasta inspect  FILE... [--framework ...] [--heavy-id ID] [--light-id ID] [--json]
```

### 8.1 seed 항체 파일을 준비하는 흐름

1. **서열을 구합니다.** 하이브리도마·파지 디스플레이로 찾은 항체의 시퀀싱
   결과(발견 플랫폼이나 CRO가 FASTA 또는 Excel로 전달), 논문·특허·데이터베이스
   (IMGT, SAbDab, UniProt, INN 목록), 또는 3.4절처럼 구조에서 뽑아낸 서열.
2. **`fasta inspect`로 어떤 레코드가 들어 있는지 봅니다.** 역할(heavy/light)이
   헤더에서 읽혔는지, 길이가 Fv 범위인지.
3. **`fasta clean`으로 정리합니다.** 줄바꿈·대소문자·갭·끝 `*`·헤더를 ABD
   형태로.
4. **`fasta validate`로 확인하고**, 같은 항체의 구조가 있다면
   `--against`로 둘이 같은 분자인지 대조합니다.

### 8.2 어느 레코드가 중쇄이고 어느 것이 경쇄인가

중쇄·경쇄를 바꿔 올리는 것이 가장 비싼 실수라서, 도구는 **서열 내용으로
추측하지 않습니다**. 역할은 다음 우선순위로만 정해집니다.

1. `--heavy-id ID` / `--light-id ID` -- 레코드 ID(헤더의 첫 단어)를 직접
   지정. 가장 확실한 방법.
2. framework가 레코드 하나를 뜻하면 그것: `nanobody` -> heavy, `scfv` -> 두
   도메인을 가진 레코드 하나.
3. 헤더의 **토큰**(글자·숫자가 아닌 문자로 쪼갠 단어). 그래서 RCSB의
   `1N8Z_2|Chains B|... heavy chain` 같은 헤더도 읽힙니다.

   | 역할 | 인식하는 토큰(대소문자 무시) |
   |---|---|
   | heavy | `H`, `HC`, `VH`, `VHH`, `heavy`, `heavychain`, `hchain`, `IGH`, `nanobody`, `Nb`, `sdAb`, 또는 `IGHV`/`IGHG`/`IGHM`/`IGHA`/`IGHE`/`IGHD`/`IGHJ`로 시작하는 유전자 이름(`IGHV3-23`) |
   | light | `L`, `LC`, `VL`, `VK`, `light`, `lightchain`, `lchain`, `kappa`, `lambda`, `IGK`, `IGL`, 또는 `IGKV`/`IGLV`/`IGKC`/`IGLC`/`IGKJ`/`IGLJ`로 시작하는 유전자 이름 |
   | scFv | `scFv`, `scFvs`, `singlechain` |

헤더에 heavy 토큰과 light 토큰이 함께 있으면 모호(`FASTA_ROLE_AMBIGUOUS`,
단 `scFv` 토큰이 있으면 scFv로 봄), 아무 토큰도 없으면
`FASTA_ROLE_UNKNOWN` -- `--framework igg`에서는 **오류**, 그 외는 경고입니다.
이때는 `--heavy-id`/`--light-id`를 주면 됩니다.

```bash
abd-prep fasta clean seq.fasta seed.fasta --framework igg --heavy-id clone7_1 --light-id clone7_2
```

### 8.3 `--framework`의 뜻

| 값 | 요구하는 레코드 | 길이 휴리스틱(경고 기준) |
|---|---|---|
| `igg` | heavy 1개 + light 1개, 역할 둘 다 확정되어야 함 | VH 105-140, VL 95-130 |
| `nanobody` | 레코드 1개 (heavy). light 토큰이 있으면 `FASTA_ROLE_CONFLICT` | VHH 100-140 |
| `scfv` | 레코드 1개 (VH + 링커 + VL) | 220-280. 링커 `(GGGGS)n` / `(GGGS)n` 위치를 정보로 표시 |

`clean`에서는 필수(출력 헤더가 역할에 따라 달라지므로), `validate`/`inspect`
에서는 선택입니다. 생략하면 레코드 수·역할 규칙은 건너뛰고
`FASTA_FRAMEWORK_UNSPECIFIED` 정보 한 줄을 남깁니다. 레코드 수로 framework를
추정하지는 않습니다.

### 8.4 `clean`이 고치는 것 -- 형식만

```bash
abd-prep fasta inspect seed_raw.fasta --framework igg
```

```
== seed_raw.fasta  records: 2
  heavy: role H (header), 120 res
     EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS
  VL_1: role L (header), 107 res
     DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK
```

`heavy`는 헤더 토큰 `heavy`로, `VL_1`은 토큰 `VL`과 `kappa`로 역할이
읽혔습니다. 이 파일은 Windows 줄바꿈, 소문자 구간, 정렬 갭 `--`, 끝 `*`를
갖고 있습니다. 정리합니다.

```bash
abd-prep fasta clean seed_raw.fasta seed.fasta --framework igg
```

```
[fasta clean] seed_raw.fasta -> seed.fasta  records: H(120 res), L(107 res)
  INFO    FASTA_GAPS_REMOVED   heavy  alignment gap character(s) '-'/'.' removed (1 run(s))
  INFO    FASTA_LOWERCASE      heavy  14 lower-case letter(s) upper-cased
  INFO    FASTA_TRAILING_STOP  VL_1   terminal '*' (translation stop) removed
  INFO    FASTA_CRLF           file   6 line(s) end with CR/LF (Windows); written as LF
  summary: 0 error(s), 0 warning(s), 4 info
  report: seed.report.json
```

출력 파일의 헤더는 ABD 규약 ID 뒤에 **원래 헤더를 설명으로 보존**합니다.

```
>H heavy trastuzumab VH, RCSB 1N8Z chain B residues 1-120 ...
EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTI...
>L VL_1 kappa trastuzumab VL, RCSB 1N8Z chain A residues 1-107 ...
DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTD...
```

`clean`이 바꾸는 것의 전부:

- CR/LF -> LF, BOM 없는 UTF-8(`--encoding`을 주면 그 인코딩에서 변환).
- 헤더를 `>H ...`, `>L ...`, `>scFv ...`로 다시 쓰고 원래 헤더를 뒤에 둠.
- 서열 안의 공백, 숫자(번호 붙은 정렬), 정렬 갭 `-`·`.` 제거, 소문자 ->
  대문자.
- **끝의** `*` 제거. **중간의** `*`는 지우지 않고 **오류**로 보고
  (`FASTA_INTERNAL_STOP`). 번역 중간에 종결 코돈이 들어간 서열이라는
  뜻이므로, 글자를 지울 일이 아니라 DNA의 읽기 틀을 확인할 일입니다.
- 모호 글자(`X`, `B`, `Z`, `J`, `U`, `O`)는 **자동으로 바꾸지 않습니다**.
  무엇으로 바꿀지 알고 있다면 `--replace-ambiguous X=G,B=N`으로 명시하고,
  그 치환은 리포트에 기록됩니다(`validate`에 주면 메모리에서만 적용되고
  파일은 `clean`만 씁니다).
- `--wrap 60`을 주면 서열을 60자마다 줄바꿈(기본은 한 줄).

그 외에는 **읽은 그대로** 씁니다. 가변 영역만 남기도록 자르지도, 잔기를
바꾸지도 않습니다.

### 8.5 `validate` 메시지 읽기

**경고 예 -- 신호 펩타이드와 불변영역이 붙어 있는 중쇄**

```
[fasta validate] seed_full.fasta
  WARNING FASTA_LENGTH_RANGE    heavy_full  207 residues is longer than a typical VH (105-140); constant region or signal peptide included? ABD seeds are variable domains
  WARNING FASTA_SIGNAL_PEPTIDE  heavy_full  starts with M followed by a hydrophobic stretch: signal peptide not removed? Mature variable domains usually start with E/Q/D
  summary: 0 error(s), 2 warning(s), 0 info
```

신호 펩타이드는 세포가 단백질을 분비할 때만 쓰고 잘라 버리는 앞부분
(보통 `M`으로 시작하는 소수성 구간 약 20잔기)입니다. 성숙한 항체 가변
도메인은 대개 `E`, `Q`, `D`로 시작합니다. 발현 벡터의 서열을 그대로 가져오면
이것이 붙어 옵니다. 도구는 **자르지 않습니다** -- 어디까지가 신호
펩타이드인지는 사용자가 판단해 출발 서열을 고쳐야 합니다.

**오류 예 -- 모호 글자와 내부 종결**

```
[fasta validate] seed_bad.fasta
  ERROR   FASTA_INTERNAL_STOP  heavy_bad  '*' inside the sequence at position(s) 61: an internal stop codon means the translation is wrong; check the DNA reading frame instead of deleting the '*'
  ERROR   FASTA_BAD_ALPHABET   heavy_bad  non-standard letter(s) X@10, B@50: only the 20 standard amino acids are accepted (X/B/Z/J/U/O are ambiguous). Replace them explicitly with --replace-ambiguous X=G or fix the source sequence
  ERROR   FASTA_BAD_ALPHABET   light_bad  non-standard letter(s) Z@20: only the 20 standard amino acids are accepted (X/B/Z/J/U/O are ambiguous). Replace them explicitly with --replace-ambiguous X=G or fix the source sequence
  summary: 3 error(s), 0 warning(s), 0 info
```

| 메시지 | 뜻 | 대처 |
|---|---|---|
| `ERROR FASTA_BAD_ALPHABET` | 표준 20종 밖의 글자(위치 표시) | 원본 서열을 고치거나, 확신이 있을 때만 `--replace-ambiguous X=G` |
| `ERROR FASTA_INTERNAL_STOP` | 서열 중간의 `*` | DNA 읽기 틀 확인. 지우지 말 것 |
| `ERROR FASTA_NUCLEIC` | A/C/G/T/U/N만으로 된 20자 이상 서열 = DNA/RNA | 단백질로 번역한 서열을 넣어야 함 |
| `ERROR FASTA_ROLE_UNKNOWN` (igg) / `WARNING` (그 외) | 헤더에서 heavy/light를 알 수 없음 | `--heavy-id`/`--light-id` |
| `ERROR FASTA_ROLE_AMBIGUOUS` | 헤더에 heavy·light 토큰이 함께 | 헤더 정리 또는 `--heavy-id`/`--light-id` |
| `ERROR FASTA_FRAMEWORK_MISMATCH` | 레코드 수·역할이 framework와 안 맞음(예: igg인데 heavy 2개) | 파일을 항체 하나씩으로 나누거나 framework 확인 |
| `WARNING FASTA_LENGTH_RANGE` | 길이가 역할의 통상 범위 밖 | 잘렸는지, 불변영역·신호 펩타이드가 붙었는지 확인 |
| `WARNING FASTA_SIGNAL_PEPTIDE` | `M` + 소수성 구간으로 시작 | 신호 펩타이드를 떼어 낸 성숙 서열인지 확인 |
| `WARNING FASTA_DUPLICATE_SEQUENCE` / `FASTA_DUPLICATE_ID` | 같은 서열·같은 ID가 둘 | 복사 실수인지 확인 |
| `INFO FASTA_SCFV_LINKER` | scFv 링커 위치 | 참고 |

휴리스틱(길이·신호 펩타이드·중복 서열)은 아주 긴 CDR-H3나 인간화 변이체에서
잘못 울릴 수 있습니다. 경고일 뿐이며 `--no-heuristics`로 끌 수 있습니다.

### 8.6 구조와 서열이 같은 항체인지 확인하기 (`--against`)

같은 항체가 두 형태로 있을 때 -- Fv 구조(`pdb split` 결과)와 서열(데이터베이스
·시퀀싱) -- ABD에 둘을 올리기 전에 서로 맞는지 확인합니다.

```bash
abd-prep fasta validate seed.fasta --framework igg --against ./out/her2-trastuzumab_fv.pdb
```

단순한 문자열 비교가 아닙니다. 결정 구조에는 보이지 않아 빠진 잔기가 있고,
seed 서열에는 태그·신호 펩타이드·불변영역 같은 여분이 붙어 있을 수 있어서,
도구는 **PDB 잔기 번호를 따라** 체인을 걸으며(번호가 뛰면 그만큼 결손,
삽입 코드는 연속으로) FASTA 안에서 가장 잘 맞는 시작 위치를 찾습니다.

| 결과 | 수준 | 뜻 |
|---|---|---|
| `FASTA_PDB_IDENTITY` | 정보 | 체인별 동일 잔기 수 / 비교한 잔기 수, 오프셋 |
| `FASTA_PDB_GAP` | 정보 | 구조에 없는(모델링되지 않은) 잔기 |
| `FASTA_EXTRA_RESIDUES` | 정보 | 서열이 구조보다 N/C 말단에서 긺(태그·불변영역) |
| `FASTA_PDB_UNCOVERED` | 경고 | 서열이 구조보다 짧음(구조에 있는 잔기가 서열에 없음) |
| `FASTA_PDB_MISMATCH` | **오류** | 둘 다 있는데 잔기가 다름 -- 점 돌연변이나 다른 클론 |
| `FASTA_PDB_NO_CORRESPONDENCE` | **오류** | 동일성 50% 미만: 다른 항체이거나 체인 대응이 틀림 |

PDB 체인이 `H`/`L`로 이름 붙어 있으면 자동으로 대응하고, 그렇지 않으면
`--chain-map`을 요구합니다. 예를 들어 원본 1N8Z(체인 A/B/C)와 비교하면:

```
abd-prep: error: cannot tell which PDB chain is heavy/light: 1N8Z.pdb has chains A, B, C and none is labelled H or L. Pass --chain-map, e.g. --chain-map H=B,L=A
```

```bash
abd-prep fasta validate seed.fasta --framework igg --against 1N8Z.pdb --chain-map H=B,L=A
```

```
[fasta validate] seed.fasta
  INFO    FASTA_PDB_IDENTITY   H:B<-H  PDB chain B vs record 'H': 120/120 aligned residues identical (PDB residue 1 = FASTA position 1)
  WARNING FASTA_PDB_UNCOVERED  H:B<-H  100 PDB residue(s) (121..220) lie outside the FASTA record; the seed sequence is shorter than the structure
  INFO    FASTA_PDB_IDENTITY   L:A<-L  PDB chain A vs record 'L': 107/107 aligned residues identical (PDB residue 1 = FASTA position 1)
  WARNING FASTA_PDB_UNCOVERED  L:A<-L  107 PDB residue(s) (108..214) lie outside the FASTA record; the seed sequence is shorter than the structure
  summary: 0 error(s), 2 warning(s), 2 info
```

원본에는 Fab 전체(불변영역 포함)가 있으므로 Fv 서열은 당연히 더 짧고,
그 사실이 경고로 나옵니다. Fv 구조(`_fv.pdb`)와 비교하면 3.4절처럼 정보
두 줄로 끝납니다. `--against`는 한 번에 FASTA 파일 하나만 받습니다.

## 9. CSV 실험 결과 표 다루기 (`abd-prep csv`)

```
abd-prep csv clean    IN.csv OUT.csv [--map STD=HEADER,...] [--delimiter CHAR] [--drop-missing] [--drop-unmapped]
abd-prep csv validate FILE... [--map STD=HEADER,...] [--delimiter CHAR] [--against seed.fasta|fv.pdb]
abd-prep csv inspect  FILE... [--map STD=HEADER,...] [--delimiter CHAR] [--json]
```

### 9.1 실험 결과 표를 준비하는 흐름

1. **Excel이라면 CSV로 저장합니다.** 파일 > 다른 이름으로 저장 > **CSV
   UTF-8(쉼표로 분리)**. 도구는 `.xlsx`를 읽지 않습니다.
2. **`csv inspect`로 열 이름과 채워진 정도를 봅니다.**
3. **`--map`으로 내 열 이름을 표준 열에 연결**하고 `csv validate`로 점검.
4. **`csv clean`으로 표준 형태(UTF-8, BOM 없음, `,` 구분, 표준 열 이름)로
   저장.**
5. seed 파일이 있으면 `--against`로 표의 행들이 그 seed를 가리키는지 확인.

### 9.2 표준 열과 `--map`

실험실마다 열 이름이 다르므로(`KD`, `Kd (nM)`, `affinity`...), `--map`으로
내 헤더를 표준 열에 투영합니다. 헤더가 이미 표준 이름과 같으면 자동으로
매핑됩니다.

| 표준 열 | 내용 | 요구 수준 |
|---|---|---|
| `id` | 클론/변이체 식별자 | **필수** (행을 구별할 수 있어야 함) |
| `sequence` | 항체 전체 서열 -- 또는 | 두 형태 중 하나를 **권장**. 없으면 `--against`는 `id`로만 대조 |
| `heavy` + `light` | 두 체인을 각각의 열에 | |
| `value` | 측정값 (KD, Tm, 발현량 ...) | 선택 |
| `unit` | `value`의 단위 (nM, pM, degC ...) | 선택. 혼용은 보고만 하고 환산하지 않음 |
| `target`, `assay`, `notes` | 항원 이름, 측정법(SPR, BLI, ELISA ...), 메모 | 선택 |

```bash
abd-prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv clean results.csv results_clean.csv --map id=Clone --map "value=KD, nM"   # 헤더에 쉼표가 있으면 따로
```

매핑되지 않은 열은 표준 열 뒤에 그대로 따라옵니다(`--drop-unmapped`로
제거). `--map`이 파일에 없는 헤더를 가리키면 실제 헤더 목록을 보여 주며
종료 2:

```
abd-prep: error: --map id=CloneID: the file has no column 'CloneID'; header is: 'Clone', 'AA_seq', 'KD', 'Unit', 'Note'
```

### 9.3 인코딩과 구분자 -- Excel이 만든 파일

**한국어 Windows Excel의 "CSV UTF-8"** 저장 결과를 `inspect`하면:

```bash
abd-prep csv inspect results.csv
```

```
== results.csv  delimiter: ;  rows: 8
  Clone: 8/8 filled  e.g. Ab001; Ab002; Ab003
  AA_seq: 7/8 filled  e.g. EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQA...; ...
  KD: 6/8 filled  e.g. 2.3; <0.1; 1,2
  Unit -> unit: 8/8 filled  e.g. nM; pM
  Note: 8/8 filled  e.g. baseline; no measurement; below assay limit
```

구분자가 `;`로 검출됐고(`clean`이 `,`로 바꿈), `Unit` 열은 표준 이름과
대소문자만 달라 자동 매핑(`Unit -> unit`)됐습니다. BOM은 `validate`에서
`FILE_BOM` 정보로 나오고 `clean`이 제거합니다.

**"CSV(쉼표로 분리)"로 저장해 CP949가 된 파일**은 UTF-8이 아니라서 그대로는
오류입니다.

```
[csv validate] results_kr.csv
  ERROR   FILE_NOT_UTF8  file  not valid utf-8: invalid start byte near byte 2; ABD requires UTF-8. If the file is in another encoding, pass it with --encoding (e.g. --encoding cp949) and run clean to convert
```

`--encoding cp949`를 주고 `clean`하면 UTF-8로 변환해 저장합니다. 열 이름이
한글이어도 `--map`에 그대로 쓸 수 있습니다.

```bash
abd-prep csv clean results_kr.csv results_kr_clean.csv --encoding cp949 --map id=클론,sequence=서열,value=KD,unit=단위
```

```
[csv clean] results_kr.csv -> results_kr_clean.csv  rows: 3
  INFO    FILE_ENCODING_CONVERTED  file  input decoded as cp949; output is written as UTF-8
  ERROR   CSV_ROW_LENGTH           rows  1 row(s) do not have the header's 5 cells: line 4 (6 cells); an unquoted delimiter inside a value is the usual cause. These rows are not written by clean
  summary: 1 error(s), 0 warning(s), 1 info
  report: results_kr_clean.report.json
```

4번째 줄의 메모에 따옴표 없는 쉼표가 들어 있어 열이 6개로 쪼개졌습니다.
이런 행은 **출력에서 빠지고** 오류로 보고되므로(종료 1), 원본에서 그 셀을
따옴표로 감싸거나 쉼표를 빼고 다시 돌리면 됩니다.

구분자 검출이 안 되면(행이 아주 적거나 따옴표가 많을 때) 헤더 줄에 후보가
하나만 있으면 그것을 쓰고, 그래도 모르면 `,`로 가정하며 경고
`CSV_DELIMITER_GUESSED`를 냅니다. `--delimiter ';'` 또는 `--delimiter tab`으로
직접 지정할 수 있습니다.

### 9.4 `validate` 메시지 읽기

```bash
abd-prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
```

```
[csv validate] results.csv
  INFO    FILE_BOM                file      UTF-8 byte-order mark present (Excel writes one); clean removes it
  INFO    CSV_DELIMITER           file      delimiter ';' (detected); clean writes ','
  WARNING CSV_MISSING_VALUES      sequence  missing values in 'sequence': 1 cell(s): line 8; --drop-missing removes these rows
  WARNING CSV_MISSING_VALUES      value     missing values in 'value': 2 cell(s): line 3, line 9; --drop-missing removes these rows
  WARNING CSV_DECIMAL_COMMA       value     1 cell(s): line 5 '1,2' contain a ',' inside the number: a European decimal mark (1,2 = 1.2) or a thousands separator (1,200)? Save the file with plain '.' decimals and no grouping -- values are not converted
  WARNING CSV_VALUE_NOT_NUMERIC   value     1 cell(s): line 4 '<0.1' are not plain numbers (assay limits like '<0.1' or text); values are kept as written
  WARNING CSV_UNIT_MIXED          unit      2 different units: 'nM' x7, 'pM' x1; values are not converted
  WARNING CSV_ID_DUPLICATE        id        1 id value(s) repeat: 'Ab001' x2
  INFO    CSV_SEQUENCE_LOWERCASE  sequence  1 cell(s): line 9 contain lower-case letters (left as-is)
  summary: 0 error(s), 6 warning(s), 3 info
```

경고는 모두 "도구가 임의로 고치지 않은 것"입니다. 값에 대한 결정은
사용자가 해야 하기 때문입니다.

| 메시지 | 뜻 | 대처 |
|---|---|---|
| `ERROR CSV_ID_MISSING` | 매핑 후 `id` 열이 없음 | `--map id=<헤더>` |
| `ERROR CSV_ROW_LENGTH` | 헤더와 셀 수가 다른 행(따옴표 없는 구분자가 흔한 원인). `clean`은 이 행을 쓰지 않음 | 원본의 그 행을 고침 |
| `ERROR CSV_NO_HEADER` / `CSV_HEADER_EMPTY_NAME` / `CSV_HEADER_DUPLICATE` | 헤더가 없거나 빈 이름·중복 이름 | 첫 줄을 고침 |
| `ERROR CSV_SEQUENCE_ALPHABET` / `CSV_SEQUENCE_NUCLEIC` | 서열 열에 표준 외 글자 / DNA처럼 보이는 값 | FASTA와 같은 규칙(8.5절) |
| `WARNING CSV_SEQUENCE_MISSING` | 서열 열이 없음 | `--against`가 `id`로만 대조됨을 알고 진행, 또는 `--map sequence=...` |
| `WARNING CSV_MISSING_VALUES` | 결측(`''`, `NA`, `N/A`, `NaN`, `None`, `null`, `-`, `--`, `ND`, `N.D.`) | 그대로 두거나 `--drop-missing`으로 해당 행 제거(표준 열의 결측만 대상) |
| `WARNING CSV_VALUE_NOT_NUMERIC` | `<0.1`, `>1000`, 문자 등 | 검출 한계 표기를 어떻게 다룰지 결정(값은 그대로 보존됨) |
| `WARNING CSV_DECIMAL_COMMA` | `1,2` -- 유럽식 소수점? 천 단위 구분? | 원본을 `.` 소수점, 자릿수 구분 없이 다시 저장 |
| `WARNING CSV_UNIT_MIXED` | 단위가 여러 가지 | 한 단위로 통일(도구는 환산하지 않음) |
| `WARNING CSV_ID_DUPLICATE` | 같은 `id`가 반복 | 반복 측정이면 의도된 것. 구분이 필요하면 id를 고침 |
| `WARNING CSV_DELIMITER_GUESSED` | 구분자 검출 실패, `,` 가정 | `--delimiter`로 지정 |

`clean`은 인코딩(UTF-8, BOM 없음), 구분자(`,`), 따옴표(최소), 헤더(표준
열을 표준 이름으로, 표준 순서로)만 바꿉니다. 셀 값은 **읽은 그대로**
씁니다 -- 숫자를 다시 쓰지도, 단위를 바꾸지도, 서열의 대소문자를 바꾸지도
않습니다.

```
id,sequence,value,unit,Note
Ab001,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPT...
Ab002,EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAAGKGLEWVARIYPT...
```

### 9.5 표의 행이 seed를 가리키는지 확인하기 (`--against`)

```bash
abd-prep csv validate results_clean.csv --against seed.fasta
```

```
  INFO    CSV_CROSSREF_MATCHED    against   3/8 row(s) match seed.fasta (0 by id, 3 by sequence)
  WARNING CSV_CROSSREF_UNMATCHED  against   5 row(s) match neither an ID nor a sequence in seed.fasta: line 3 id 'Ab002'; line 4 id 'Ab003'; line 5 id 'Ab004'; line 6 id 'Ab005'; line 8 id 'Ab006'
```

행은 `id`가 레코드/체인 ID와 같거나, `sequence`(또는 `heavy`+`light` 둘 다)가
레코드/체인 서열과 같으면 매칭입니다. 참조는 seed FASTA 또는 Fv PDB 둘 다
됩니다. 위 예에서는 변이체 5개가 seed와 서열이 달라 매칭되지 않았는데,
변이체 표라면 당연한 결과입니다 -- 이 검사는 "seed 자체가 표에 들어 있나",
"id 체계가 서로 맞나"를 보는 용도입니다.

## 10. 스크립트와 CI에서 쓰기

종료 코드와 finding code가 고정되어 있으므로 셸 스크립트나 CI(지속적 통합,
자동 검사 파이프라인)에서 쓸 수 있습니다.

```bash
#!/usr/bin/env bash
set -euo pipefail

# 1. 구조 두 파일 생성. 경고도 실패로 보려면 --strict (1N8Z는 해상도 경고 때문에 실패함)
abd-prep pdb split --preset her2-trastuzumab --out ./out

# 2. seed 서열 정리 + 구조와 대조. 오류가 있으면 종료 코드 1로 스크립트가 멈춘다
abd-prep fasta clean seed_raw.fasta ./out/seed.fasta --framework igg
abd-prep fasta validate ./out/seed.fasta --framework igg --against ./out/her2-trastuzumab_fv.pdb

# 3. 결과 표. 경고 목록만 뽑아 보기 (jq가 있을 때)
abd-prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit --json \
  | jq -r '.findings[] | select(.level=="warning") | "\(.code)\t\(.where)\t\(.message)"'
```

- `--json`은 리포트 전체(6.4절의 구조)를 표준 출력으로 냅니다. `findings[]`의
  `code`로 걸러 내세요.
- `--strict`는 모든 경고를 실패로 바꿉니다. 1N8Z의 `PDB_RESOLUTION`처럼
  사실이지만 어쩔 수 없는 경고도 포함되므로, 기준 파일에서 한 번 돌려 보고
  어떤 경고가 나오는지 확인한 뒤 적용하세요.
- `validate`는 파일 여러 개를 한 번에 받고, 메시지 위치 앞에 파일 이름을
  붙여 구분합니다. 종료 코드는 가장 나쁜 파일의 것입니다.
- `--report PATH`로 리포트를 원하는 위치에 두고, 결과 파일과 함께 보관하세요.

Python에서 바로 부르려면 `abd_prep.cli.main([...])`에 인자 목록을 넘기면
종료 코드가 반환됩니다. 더 깊은 사용(라이브러리로 쓰기)은 개발자 가이드를
보세요.

## 11. 제한 사항

- **구형 PDB 형식만** 읽습니다(고정 칸). mmCIF는 읽지 않습니다. RCSB가
  mmCIF로만 제공하는 아주 큰 항목은 `fetch`가 그렇다고 알려 줍니다.
- 출력에는 `ATOM`, `HETATM`, `TER`, `END`만 씁니다. `HEADER`, `REMARK`,
  `SEQRES`, `SSBOND`, `LINK`, `CRYST1`, `HELIX`/`SHEET` 등은 버립니다.
  출처(항목 ID, 개정일, SHA-256)는 대신 리포트에 남습니다.
- 대체 위치는 **원자 단위로 점유율 기준** 선택입니다. 잔기 전체를 한
  컨포머로 맞추지 않으며, PyMOL의 "A 남기기"와 결과가 다를 수 있습니다.
- **Fv 경계를 계산하지 않습니다.** 프리셋은 트라스투주맙의 문헌값이고, 다른
  항체는 범위를 직접 지정해야 합니다. FASTA 그룹은 신호 펩타이드·불변영역이
  붙었을 가능성을 경고하지만 자르지 않습니다. ANARCI 같은 항체 번호 매기기
  도구는 포함되어 있지 않습니다.
- FASTA 휴리스틱(길이 범위, 신호 펩타이드)은 특이한 항체(매우 긴 CDR-H3,
  인간화 변이체)에서 잘못 울릴 수 있습니다. 경고이며 `--no-heuristics`로
  끌 수 있습니다.
- `--against`는 체인당 오프셋 하나로 맞춥니다. "같은 분자, 잔기 몇 개
  결손" 용도이고, 서로 다른 항체를 정렬하는 도구가 아닙니다(그 경우 "대응
  없음" 한 줄로 끝남).
- CSV: `.xlsx`를 읽지 않습니다(Excel에서 CSV UTF-8로 저장). 값은 **형식**만
  검사합니다 -- 단위 환산, 범위·타당성 검사는 없습니다. 셀 안의 줄바꿈은
  처리하지만 헤더 줄은 반드시 있어야 합니다.
- `--strict`는 모든 경고를 실패로 바꿉니다(위의 휴리스틱과 `PDB_RESOLUTION`
  포함).
- **ABD의 검증을 대신하지 않습니다.** 종료 코드 0은 "알려진 형식 문제가
  없다"는 뜻이고, ABD가 그 생물학적 내용을 받아 준다는 보장은 아닙니다.

## 12. 아직 확인되지 않은 ABD 규약 3가지

공개된 ABD 사용자 가이드는 요구 사항을 "proper chain identification
(H/L/T)", "verify sequence format and headers" 정도로만 적고 있습니다. 그래서
세 가지는 **가장 그럴듯한 해석**으로 구현되어 있고, 각각 파일 한 곳에
격리되어 있어 확인되면 바로 고칠 수 있습니다.

| 규약 | 현재 동작 | 만약 ABD가 다르게 요구하면 |
|---|---|---|
| PDB 타겟·항체 파일의 체인 ID | Fv는 `L`/`H`로 쓰고, 타겟은 원본 ID 유지 | `--rename-target T`(또는 `pdb clean --rename C=T`)로 즉시 대응 가능 |
| FASTA 헤더 | `>H ...`, `>L ...`, `>scFv ...` + 원래 헤더를 설명으로 | 헤더가 거부되면 `>H`/`>L`만 남기고 설명을 지우면 됨. 도구 쪽 수정은 상수 한 곳(개발자 가이드 참조) |
| CSV 열 스키마 | 9.2절의 표준 열(도구가 정한 것) | ABD가 요구하는 열 이름으로 `--map`을 거꾸로 쓰거나 헤더만 바꾸면 됨 |

콘솔에서 이 중 하나를 확인하셨다면 저장소에 이슈로 알려 주세요.

## 13. 문제 해결 (FAQ)

**Q. `abd-prep: command not found`**
설치하지 않은 상태입니다. 저장소 디렉터리에서 `python3 -m abd_prep ...`로
실행하거나 2.2절처럼 설치하세요.

**Q. `No module named abd_prep`**
`python3 -m abd_prep`는 저장소 디렉터리(또는 설치된 환경) 안에서만 됩니다.
다른 디렉터리에서 쓰려면 설치하거나 `PYTHONPATH=/path/to/abd-prep`을
붙이세요.

**Q. "`results.xlsx` is an Excel workbook, which this tool cannot read"**
Excel에서 파일 > 다른 이름으로 저장 > CSV UTF-8(쉼표로 분리)로 저장한
뒤 그 `.csv`를 넣으세요.

**Q. `FILE_NOT_UTF8 ... invalid start byte near byte 2`**
파일이 UTF-8이 아닙니다(한국어 Windows면 대개 CP949). `--encoding cp949`를
붙여 `clean`하면 UTF-8로 변환됩니다. 어떤 인코딩인지 모르면
`--encoding latin-1`로 일단 읽어 글자가 깨지는지 보세요.

**Q. `PDB_CHAIN_NOT_FOUND ... chain 'Z' is not in the file (chains present: ABC)`**
체인 ID를 잘못 썼습니다. 메시지에 실제 체인 목록이 나옵니다.
`pdb inspect`로 먼저 확인하세요.

**Q. `PDB_SELECTION_EMPTY`**
범위가 아무 원자도 고르지 않았습니다. 잔기 번호가 파일에 적힌 번호인지
(순차 번호가 아님) 확인하세요. 범위 안에 HETATM만 있어서 전부 지워진
경우에는 `--keep-het` 힌트가 함께 나옵니다.

**Q. `range '10-5': start is after end`**
`A:5-10`처럼 작은 번호를 먼저 쓰세요. 종료 코드 2(옵션 오류)입니다.

**Q. `clean`이 오류를 냈는데 출력 파일이 생겼습니다. 써도 되나요?**
결과를 눈으로 보라고 쓴 것입니다. 종료 코드 1인 파일을 ABD에 올리지는
마세요. 오류를 고친 뒤 다시 돌리세요.

**Q. `FASTA_ROLE_UNKNOWN`이 뜹니다. 헤더를 어떻게 써야 하나요?**
가장 쉬운 방법은 `--heavy-id`/`--light-id`로 레코드 ID를 지정하는 것입니다.
헤더를 고치려면 `>clone7 heavy`처럼 8.2절의 토큰 하나를 넣으면 됩니다.

**Q. 1N8Z 프리셋이 `--strict`에서 실패합니다.**
해상도 2.52 Å 경고(`PDB_RESOLUTION`) 때문입니다. 사실이므로 도구는 숨기지
않습니다. `--strict` 없이 쓰거나, 더 좋은 해상도의 구조를 고르세요.

**Q. `--against`가 `FASTA_PDB_NO_CORRESPONDENCE`를 냅니다.**
두 파일이 다른 항체이거나 `--chain-map`의 대응(어느 PDB 체인이 heavy인지)이
뒤바뀐 경우입니다. `pdb inspect`로 각 체인의 서열을 보고 FASTA와 눈으로
맞춰 보세요.

**Q. 프록시 환경에서 `fetch`가 안 됩니다.**
`HTTPS_PROXY` 환경 변수를 설정하세요. 네트워크가 막힌 곳이면 다른
컴퓨터에서 받은 `.pdb` 파일을 `--input`으로 넣으면 `fetch` 없이 됩니다.

**Q. report.json에 내 홈 디렉터리 경로가 들어가나요?**
기본값은 파일 이름만 적습니다(`--report-paths basename`). 절대 경로가
필요하면 `--report-paths full`.

## 14. 용어집

| 용어 | 뜻 |
|---|---|
| **ABD (Amazon Bio Discovery)** | AWS의 항체 설계 서비스. 이 도구는 그 입력 파일을 준비하는 독립 도구이며 ABD의 구성 요소가 아님 |
| **레시피(recipe)** | ABD의 작업 묶음. De novo 생성, 유도 진화, 속성 예측 등 |
| **seed 항체** | 개량·평가의 출발점이 되는 기존 항체(서열 또는 구조) |
| **타겟 항원(target antigen)** | 항체가 붙어야 할 목표 단백질. ABD에는 PDB 구조로 올림 |
| **De novo 설계** | seed 없이 항원 구조만 보고 새 항체를 생성 |
| **Directed evolution(유도 진화)** | 기존 항체를 조금씩 바꿔 가며 개량 |
| **Property prediction(속성 예측)** | 결합력·안정성·개발 적합성을 계산으로 예측 |
| **Fine-tuning(미세조정)** | 자기 실험 데이터로 모델을 자기 과제에 맞게 추가 학습 |
| **골격(framework)** | 항체의 형식. IgG(중쇄+경쇄), Nanobody/VHH(단일 도메인), scFv(두 도메인을 링커로 이은 단일 사슬) |
| **중쇄(heavy, H) / 경쇄(light, L)** | 항체를 이루는 두 종류의 폴리펩타이드 사슬 |
| **Fab / Fv / VH / VL** | Fab = 항원 결합 조각(가변+불변 일부). Fv = 가변 영역 두 도메인(VH + VL)만. ABD 항체 파일은 보통 Fv |
| **불변영역(constant region)** | 항체에서 항원 결합과 무관한 뒷부분. seed에는 보통 넣지 않음 |
| **신호 펩타이드(signal peptide)** | 분비 과정에서 잘려 나가는 N-말단 앞부분(약 20잔기, `M`으로 시작). 성숙 서열에는 없어야 함 |
| **CDR** | 항원과 직접 닿는 가변 영역 안의 고리 6개. CDR-H3가 가장 다양함 |
| **링커(linker)** | scFv에서 VH와 VL을 잇는 유연한 짧은 서열(보통 `GGGGS` 반복) |
| **Kabat 번호** | 항체 잔기 번호 체계. 삽입 코드(`82A`)를 써서 길이가 다른 항체를 같은 번호로 맞춤 |
| **PDB (Protein Data Bank / 파일 형식)** | 단백질 구조 데이터베이스이자 그 고정폭 텍스트 형식의 이름 |
| **RCSB** | PDB를 운영하는 미국 기관(rcsb.org). 공개 구조를 내려받는 곳 |
| **mmCIF** | PDB 형식의 후속 형식. 이 도구는 읽지 않음 |
| **ATOM / HETATM** | 단백질 원자 / 그 외(물·이온·당·리간드·변형 잔기) 레코드 |
| **체인(chain)** | 한 폴리펩타이드 사슬. 글자 하나(`A`, `L`, `H`)로 구분 |
| **잔기(residue)** | 사슬을 이루는 아미노산 하나. 세 글자 이름 + 번호 (+ 삽입 코드) |
| **대체 위치(altloc)** | 한 원자에 좌표가 두 세트 이상 있는 상태. ABD가 거부함 |
| **점유율(occupancy)** | 대체 위치 각각이 실제로 그 자리에 있을 비율(합 1) |
| **컨포머(conformer)** | 같은 분자의 다른 형태(여기서는 대체 위치의 각 좌표 세트) |
| **TER / END** | 체인 끝 표시 / 파일 끝 표시 |
| **MODEL** | 모델이 여러 개인 파일(NMR 등)에서 각 모델의 묶음 |
| **해상도(resolution, Å)** | 구조의 선명도. 작을수록 좋음. 1 Å = 0.1 nm |
| **CA(알파 탄소)** | 각 아미노산의 중심 탄소. 연속 잔기의 CA 사이는 약 3.8 Å |
| **변형 잔기(modified residue)** | MSE, CSO 등 화학적으로 바뀐 아미노산. HETATM으로 적히지만 사슬의 일부 |
| **글리칸(glycan), NAG** | 단백질에 붙은 당. NAG는 그 구성 단위 중 하나 |
| **FASTA** | `>헤더` + 서열 줄로 된 서열 텍스트 형식 |
| **레코드(record) / ID** | FASTA의 헤더+서열 한 묶음 / 헤더의 첫 단어 |
| **1문자 기호** | 아미노산 20종을 한 글자로(A, C, D, ..., Y) |
| **모호 글자(ambiguous)** | `X`(모름), `B`(D 또는 N), `Z`(E 또는 Q), `J`, `U`, `O`. ABD가 거부 |
| **종결 코돈(`*`)** | DNA 번역이 끝나는 자리 표시. 끝에 있으면 무해, 중간에 있으면 서열 오류 |
| **CSV** | 쉼표(또는 `;`, 탭)로 열을 나눈 표 텍스트 |
| **구분자(delimiter)** | 열을 나누는 문자 |
| **헤더(CSV)** | 첫 줄의 열 이름 |
| **UTF-8 / CP949(EUC-KR)** | 글자를 바이트로 적는 방식. ABD는 UTF-8만. 한국어 Windows Excel의 기본은 CP949 |
| **BOM** | UTF-8 파일 앞의 보이지 않는 3바이트 표식. Excel이 붙임 |
| **결측(missing value)** | 값이 없음을 나타내는 셀(`NA`, `N/A`, `-` 등) |
| **KD** | 해리 상수. 결합 세기. 작을수록 강하게 붙음(nM, pM) |
| **Tm** | 녹는 온도. 열 안정성 지표(degC) |
| **SPR / BLI / ELISA** | 결합을 측정하는 실험 기법들 |
| **CRO** | 수탁 연구 기관. 실험을 대행하고 결과를 Excel 등으로 전달 |
| **finding / finding code** | 도구가 보고하는 발견 사항 / 그 고정 식별자(`PDB_ALTLOC_PRESENT` 등) |
| **종료 코드(exit code)** | 명령이 끝날 때 운영체제에 돌려주는 숫자. 0 통과, 1 검증 오류, 2 입력 오류 |
| **SHA-256** | 파일 내용의 지문. 한 바이트만 달라도 값이 바뀜 |
| **CI** | 코드나 데이터가 바뀔 때마다 자동으로 검사하는 파이프라인 |

## 부록 A. 명령 치트시트

```bash
# ---------- 한 번에 (간단 모드) ----------
abd-prep quickstart                                           # ./abd-quickstart/ 에 타겟 PDB·Fv PDB·seed FASTA·CSV 예시·README
abd-prep quickstart --out ./my-test --preset her2-trastuzumab-6lbx-6b9z --force

# ---------- PDB ----------
abd-prep pdb fetch 1N8Z 6LBX --out ./raw                     # RCSB에서 받기(캐시 ~/.cache/abd-prep)
abd-prep pdb split --preset her2-trastuzumab --out ./out      # 프리셋으로 타겟+Fv 두 파일
abd-prep pdb split --pdb 1N8Z --target C:489-607 --heavy B:1-120 --light A:1-107 \
                   --rename-heavy H --rename-light L --out ./out
abd-prep pdb split --input my.pdb --target B:530-601 --name her2_frag --out ./out
abd-prep pdb clean in.pdb out.pdb --keep A:1-107 --keep B:1-120 --rename A=L,B=H
abd-prep pdb clean in.pdb out.pdb --keep B --std-residues --strip-h --renumber
abd-prep pdb validate out.pdb --expect-chains 2
abd-prep pdb validate out_target.pdb out_fv.pdb --preset her2-trastuzumab
abd-prep pdb inspect out.pdb
abd-prep pdb inspect out.pdb --fasta > seed.fasta

# ---------- FASTA ----------
abd-prep fasta inspect seed_raw.fasta --framework igg
abd-prep fasta clean seed_raw.fasta seed.fasta --framework igg
abd-prep fasta clean seq.fasta seed.fasta --framework igg --heavy-id clone7_1 --light-id clone7_2
abd-prep fasta clean vhh.fasta seed.fasta --framework nanobody
abd-prep fasta clean scfv.fasta seed.fasta --framework scfv --wrap 60
abd-prep fasta validate seed.fasta --framework igg
abd-prep fasta validate seed.fasta --framework igg --against fv.pdb
abd-prep fasta validate seed.fasta --framework igg --against 1N8Z.pdb --chain-map H=B,L=A
abd-prep fasta validate seed.fasta --framework igg --replace-ambiguous X=G --no-heuristics

# ---------- CSV ----------
abd-prep csv inspect results.csv
abd-prep csv validate results.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv clean results.csv results_clean.csv --map id=Clone,sequence=AA_seq,value=KD,unit=Unit
abd-prep csv clean results_kr.csv results_clean.csv --encoding cp949 --map id=클론,sequence=서열,value=KD,unit=단위
abd-prep csv clean results.csv results_clean.csv --map id=Clone --delimiter ';' --drop-missing --drop-unmapped
abd-prep csv validate results_clean.csv --against seed.fasta

# ---------- 공통 ----------
abd-prep <group> <action> ... --json            # 리포트를 JSON으로
abd-prep <group> <action> ... --strict          # 경고도 실패로
abd-prep <group> <action> ... --report r.json   # 리포트 위치 지정
abd-prep <group> <action> ... --no-report       # 리포트 안 씀
abd-prep <group> <action> --help                # 전체 옵션
```

## 부록 B. 자주 만나는 finding code 빠른 대처표

전체 목록(98개)은 README의 [Finding codes](../README.md#finding-codes)에
있습니다. 아래는 실제로 자주 만나는 것만 골랐습니다.

| 코드 | 수준 | 한 줄 대처 |
|---|---|---|
| `FILE_TOO_LARGE` | 오류 | 100 MB 미만으로. 구조는 체인 선택, 표는 행 분할 |
| `FILE_NOT_UTF8` | 오류 | `--encoding cp949` 등으로 `clean` |
| `FILE_BOM` | 정보 | `clean`이 제거 |
| `FILE_EXTENSION` | 경고 | 확장자가 어색하지만 파싱은 함. `.xlsx`만 입력 오류 |
| `PDB_ALTLOC_PRESENT` | 오류 | `pdb clean` |
| `PDB_HETATM_PRESENT` | 경고 | `pdb clean` (의도적이면 `--allow-het`) |
| `PDB_TER_*` | 오류 | `pdb clean`이 재생성 |
| `PDB_GAP_CA_DISTANCE` | 경고 | 결정에 없는 고리. 사실이면 그대로 |
| `PDB_RESIDUE_COUNT` | 경고 | Fv 범위 재확인 |
| `PDB_MODIFIED_RESIDUE_REMOVED` | 경고 | `clean --std-residues` |
| `PDB_RESOLUTION` | 경고 | 사실 보고. 다른 구조 검토 |
| `PDB_CHAIN_NOT_FOUND` / `PDB_SELECTION_EMPTY` | 오류 | `pdb inspect`로 체인·번호 확인 |
| `FASTA_ROLE_UNKNOWN` / `FASTA_ROLE_AMBIGUOUS` | 오류·경고 | `--heavy-id`/`--light-id` |
| `FASTA_BAD_ALPHABET` | 오류 | 원본 수정 또는 명시적 `--replace-ambiguous` |
| `FASTA_INTERNAL_STOP` | 오류 | DNA 읽기 틀 확인. 지우지 말 것 |
| `FASTA_NUCLEIC` | 오류 | 단백질 서열로 번역해서 넣기 |
| `FASTA_FRAMEWORK_MISMATCH` | 오류 | 레코드 수·framework 확인 |
| `FASTA_LENGTH_RANGE` / `FASTA_SIGNAL_PEPTIDE` | 경고 | 신호 펩타이드·불변영역 포함 여부 확인 |
| `FASTA_PDB_MISMATCH` | 오류 | 서열과 구조가 다른 잔기. 어느 쪽이 맞는지 확인 |
| `FASTA_PDB_NO_CORRESPONDENCE` | 오류 | 다른 항체이거나 `--chain-map` 뒤바뀜 |
| `CSV_ID_MISSING` | 오류 | `--map id=<헤더>` |
| `CSV_ROW_LENGTH` | 오류 | 해당 행의 따옴표 없는 구분자 수정 |
| `CSV_MISSING_VALUES` | 경고 | 그대로 두거나 `--drop-missing` |
| `CSV_VALUE_NOT_NUMERIC` / `CSV_DECIMAL_COMMA` | 경고 | 값 표기 정리(`<0.1`, `1,2`) |
| `CSV_UNIT_MIXED` | 경고 | 단위 통일(도구는 환산 안 함) |
| `CSV_CROSSREF_UNMATCHED` | 경고 | id 체계·서열이 seed와 맞는지 확인 |
