# abd-prep 쉬운 사용 방법 (Easy Tutorial)

> **20분 안에** Amazon Bio Discovery(ABD) 테스트에 필요한 파일을 전부 만들어
> 봅니다. 터미널(명령을 글자로 입력하는 검은 창)을 거의 써 본 적이 없어도
> 따라올 수 있게 썼습니다. 한 번에 다 되는 `quickstart` 명령 하나가 핵심이고,
> 나머지는 "그 다음에 무엇을 하면 되는지"입니다.
>
> 더 자세한 설명은 [사용자 가이드](user-guide.ko.md), 명령·옵션의 기준 명세는
> 영어 [README](../README.md)에 있습니다.

## 목차

- [1. 무엇을 하는 도구인가](#1-무엇을-하는-도구인가)
- [2. 준비물 확인 (2분)](#2-준비물-확인-2분)
- [3. 내려받기 (1분)](#3-내려받기-1분)
- [4. 한 번에 준비하기: `quickstart` (1분)](#4-한-번에-준비하기-quickstart-1분)
- [5. ABD에 올리기: 어떤 파일을 어디에](#5-abd에-올리기-어떤-파일을-어디에)
- [6. 내 파일로 바꾸기](#6-내-파일로-바꾸기)
- [7. 메시지 읽는 법](#7-메시지-읽는-법)
- [8. 자주 막히는 곳](#8-자주-막히는-곳)
- [9. 다음에 읽을 것](#9-다음에-읽을-것)
- [용어 열 개](#용어-열-개)

## 1. 무엇을 하는 도구인가

ABD에 항체 설계 작업을 맡기려면 파일을 올려야 합니다. **타겟 항원**(항체가
붙어야 하는 표적 단백질)의 **구조 파일**(PDB), 출발점이 되는 **seed 항체**의
구조 또는 **서열 파일**(FASTA, 아미노산을 글자 한 개씩 적은 텍스트), 그리고
실험 결과를 담은 **표 파일**(CSV)입니다. 이 파일들은 뷰어에서는 멀쩡해
보여도 ABD의 파서(파일을 읽어 들이는 프로그램)는 거부하는 경우가 많습니다.
한 원자가 두 자리에 기록된 "대체 위치", 지워진 체인의 `TER` 줄, 중쇄·경쇄
구분이 없는 서열 헤더, Excel이 붙여 넣는 눈에 안 보이는 글자 같은 것들입니다.

`abd-prep`은 그 **형식 문제만** 고치고, 나머지는 전부 메시지로 알려 줍니다.
좌표·서열·측정값 같은 **과학적 내용은 절대 바꾸지 않습니다.**

## 2. 준비물 확인 (2분)

필요한 것은 **Python 3.8 이상**과 인터넷 연결(구조 파일을 RCSB 공개
데이터베이스에서 내려받기 위해)입니다. 다른 설치는 없습니다.

**터미널 열기**

| 운영체제 | 방법 |
|---|---|
| macOS | `Cmd + Space` -> "터미널" 입력 -> Enter |
| Windows | 시작 메뉴 -> "PowerShell" 입력 -> Enter |
| Linux | `Ctrl + Alt + T` |

**Python 확인** -- 아래를 입력하고 Enter.

```bash
python3 --version
```

`Python 3.9.25`처럼 3.8 이상이 나오면 됩니다. Windows에서 `python3`이
없다고 하면 `py -3 --version`을 써 보고, 이후 모든 명령에서 `python3`을
`py -3`으로 바꿔 읽으세요. 아무것도 없으면
[python.org](https://www.python.org/downloads/)에서 설치합니다(설치 화면의
"Add Python to PATH"에 체크).

## 3. 내려받기 (1분)

`git`이 있으면:

```bash
git clone https://github.com/hmkim/abd-prep.git
cd abd-prep
```

`git`이 없으면 저장소 페이지의 **Code -> Download ZIP**으로 받아 풀고,
터미널에서 그 폴더로 들어갑니다(`cd` 뒤에 폴더 경로를 적습니다).

설치 과정은 없습니다. 폴더 안에서 바로 실행합니다.

```bash
python3 -m abd_prep --version
```

`abd-prep 0.1.0`이 나오면 준비 끝입니다.

## 4. 한 번에 준비하기: `quickstart` (1분)

```bash
python3 -m abd_prep quickstart
```

이 한 줄이 하는 일: RCSB에서 1N8Z(HER2 항원과 트라스투주맙 항체가 결합한
결정 구조)를 내려받고, 항원 부분과 항체 Fv(항원을 붙잡는 가변 영역) 부분을
잘라 각각 정리·검증하고, 항체 서열을 FASTA로 뽑아 구조와 대조하고, 결과
표의 예시를 만들고, 각 파일이 무엇인지 적은 README를 씁니다. 실제 출력은
이렇습니다.

```
[quickstart] her2-trastuzumab  ->  abd-quickstart
  WARNING PDB_RESOLUTION        1N8Z                       resolution 2.52 A is coarser than the 2.5 A the ABD guide recommends
  INFO    PDB_HETATM_REMOVED    target/file                removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    PDB_GAP_NUMBERING     target/C:580-591           residue numbering jumps 580 -> 591 (expected)
  INFO    PDB_GAP_CA_DISTANCE   target/C:580-591           CA-CA distance 6.7 A between consecutive residues (> 4.2 A): chain break -- expected for this preset
  INFO    PDB_HETATM_REMOVED    fv/file                    removed HETATM records: HOH x79, NAG x28, SO4 x5
  INFO    FASTA_PDB_IDENTITY    seed/H:H<-antibody_fv.pdb  PDB chain H vs record 'antibody_fv.pdb': 120/120 aligned residues identical (PDB residue 1 = FASTA position 1)
  INFO    FASTA_PDB_IDENTITY    seed/L:L<-antibody_fv.pdb  PDB chain L vs record 'antibody_fv.pdb': 107/107 aligned residues identical (PDB residue 1 = FASTA position 1)
  INFO    CSV_CROSSREF_MATCHED  csv/against                1/1 row(s) match antibody_seed.fasta (0 by id, 1 by sequence)
  summary: 0 error(s), 1 warning(s), 7 info
  report: abd-quickstart/quickstart.report.json

Files in abd-quickstart:
  target.pdb           antigen structure (chain C 109 res)
                       -> target antigen -- de novo design (required); directed evolution with antigen (optional)
  antibody_fv.pdb      seed antibody structure (chain L 107 res, chain H 120 res)
                       -> seed antibody as a structure -- directed evolution, property prediction
  antibody_seed.fasta  seed antibody sequence (>H 120 res, >L 107 res)
                       -> seed antibody as a sequence -- directed evolution, property prediction (same antibody as antibody_fv.pdb)
  assay_example.csv    assay-table template (1 placeholder row, not a measurement)
                       -> assay / results table template -- fine-tuning, results analysis; replace the placeholder row with your data
  README.txt           what each file is, where it came from, what to do next
Next: read README.txt, then create an experiment in the ABD console and upload the file(s) its recipe asks for.
```

읽는 요령 세 가지면 충분합니다.

- **`summary: 0 error(s)`** -- 오류가 0이면 성공입니다. 경고(WARNING) 1건은
  "1N8Z의 해상도 2.52 Å가 ABD 권장치 2.5 Å보다 아주 조금 거칠다"는 사실
  보고이고, 파일을 못 쓰게 하지는 않습니다.
- **INFO** 줄은 "이런 일을 했다"는 기록입니다. 물(HOH)·당(NAG)·황산 이온(SO4)을
  지웠고, 타겟 체인 C의 581-590번 고리는 원래 결정에서 보이지 않아 비어 있으며
  (`expected`), FASTA와 구조의 서열이 120/120, 107/107로 똑같다는 내용입니다.
- **`Files in abd-quickstart:`** 아래가 만들어진 파일과 쓰는 곳입니다.

새 폴더 `abd-quickstart/`에는 여섯 파일이 있습니다.

| 파일 | 무엇인가 | 어디에 쓰나 |
|---|---|---|
| `target.pdb` | 타겟 항원(HER2 도메인 IV) 구조, 체인 C 109잔기 | De novo 설계의 **타겟 항원**; 항원을 포함하는 Directed evolution |
| `antibody_fv.pdb` | seed 항체(트라스투주맙 Fv) 구조, 체인 L·H | Directed evolution / Property prediction의 **seed 항체**(구조로 올릴 때) |
| `antibody_seed.fasta` | 같은 항체의 서열(`>H` 120잔기, `>L` 107잔기) | 같은 두 레시피의 **seed 항체**(서열로 올릴 때) |
| `assay_example.csv` | 결과 표의 **형식 예시** -- 한 줄짜리, 값은 자리 표시용 | 미세조정(fine-tuning)·결과 분석에 올릴 표의 모양 참고. **그대로 올리는 파일이 아닙니다** |
| `README.txt` | 위 내용 + 출처(RCSB 개정일) + 다음 단계 | 파일과 같이 보관 |
| `quickstart.report.json` | 무엇을 했는지, 각 파일의 SHA-256(내용 지문) | 재현·추적용. 파일과 같이 보관 |

같은 명령을 다시 실행하면 "이미 파일이 있다"며 멈춥니다(실수로 덮어쓰지
않도록). 다시 만들려면 `--force`를, 다른 폴더에 만들려면 `--out 폴더이름`을
붙입니다.

```bash
python3 -m abd_prep quickstart --out ./my-first-abd-test
python3 -m abd_prep quickstart --force
```

## 5. ABD에 올리기: 어떤 파일을 어디에

ABD는 작업을 **레시피**(recipe)로 고릅니다. 하고 싶은 일에서 출발해
찾으세요.

| 하고 싶은 일 | ABD 레시피 | 올릴 파일 |
|---|---|---|
| 이 항원에 붙는 **새 항체를 처음부터** 만들어 보기 | Antibody generation -> **De novo design** | `target.pdb` (타겟 항원). 골격(Nanobody / IgG / scFv)은 콘솔에서 선택 |
| 기존 항체(트라스투주맙)를 **개량**해 보기 | Antibody generation -> **Directed evolution** | `antibody_seed.fasta` 또는 `antibody_fv.pdb` (seed). 항원을 함께 넣는 모드면 `target.pdb`도 |
| 기존 항체의 **결합력·안정성을 예측**해 보기 | **Property prediction** | `antibody_seed.fasta` 또는 `antibody_fv.pdb` (seed) |
| 내 실험 결과로 **미세조정**하거나 결과를 **대조** | fine-tuning / results | `assay_example.csv`의 열 구조를 본 뒤, **내 측정값**으로 채운 CSV |

올리기 전에 세 가지만 기억하세요.

1. `assay_example.csv`의 유일한 행은 seed 항체 자신과 자리 표시용 숫자
   `1.0 nM`입니다. 측정값이 아닙니다(`notes` 칸에도 그렇게 적혀 있습니다).
   내 데이터로 바꾼 뒤에 올립니다(6장).
2. 타겟 파일의 체인 ID는 원본 그대로 `C`입니다. 콘솔이 특정 ID(예: `T`)를
   요구하면 한 줄로 바꿀 수 있습니다.
   ```bash
   python3 -m abd_prep pdb clean abd-quickstart/target.pdb abd-quickstart/target_T.pdb --rename C=T
   ```
3. 콘솔이 파일을 거부하면 메시지를 그대로 들고 사용자 가이드
   [12장(아직 확인되지 않은 규약 3가지)](user-guide.ko.md#12-아직-확인되지-않은-abd-규약-3가지)를
   보세요. 체인 ID, FASTA 헤더, CSV 열 이름은 공개 문서에 글자 단위 명세가
   없어 "가장 그럴듯한 해석"으로 돼 있고, 각각 한 줄로 고칠 수 있습니다.

## 6. 내 파일로 바꾸기

예제가 돌아가면 다음은 내 데이터입니다. 상황별로 명령 한두 줄입니다. 모든
명령은 입력 파일을 바꾸지 않고 새 파일을 씁니다.

**내 구조 파일(PDB)이 있다.** 어떤 체인이 있는지 먼저 보고, 필요한 체인만
골라 정리합니다.

```bash
python3 -m abd_prep pdb inspect my_complex.pdb
python3 -m abd_prep pdb clean my_complex.pdb my_fv.pdb --keep A:1-107 --keep B:1-120 --rename A=L,B=H
python3 -m abd_prep pdb clean my_complex.pdb my_target.pdb --keep C
```

`--keep A:1-107`은 "체인 A의 1번부터 107번 잔기까지만", `--rename A=L,B=H`는
"체인 이름을 A는 L로, B는 H로"라는 뜻입니다. 범위를 모르면 `inspect`가
보여 주는 체인별 첫·끝 번호를 참고하세요.

**내 항체 서열(FASTA)이 있다.** 논문·특허·데이터베이스에서 복사한 서열이
흔합니다. 어느 레코드가 중쇄(heavy)이고 경쇄(light)인지 헤더에 적혀 있지
않으면 직접 알려 줍니다.

```bash
python3 -m abd_prep fasta clean paper.fasta seed.fasta --framework igg
python3 -m abd_prep fasta clean paper.fasta seed.fasta --framework igg --heavy-id VH_clone7 --light-id VL_clone7
```

`--framework igg`는 "중쇄 1개 + 경쇄 1개짜리 보통 항체"라는 뜻입니다.
나노바디(중쇄 하나짜리 단일 도메인 항체)는 `nanobody`, 두 도메인을 링커로
이어 붙인 단일 사슬이면 `scfv`입니다.

**내 실험 결과 표(CSV/Excel)가 있다.** Excel에서는 **"다른 이름으로 저장 ->
CSV UTF-8"**로 저장하고, 내 열 이름이 표준 열(id, sequence 또는
heavy+light, value, unit ...)의 무엇에 해당하는지 `--map`으로 알려 줍니다.

```bash
python3 -m abd_prep csv clean results.csv results_clean.csv --map id=Clone,heavy=VH_seq,light=VL_seq,value=KD,unit=Unit
python3 -m abd_prep csv validate results_clean.csv --against seed.fasta
```

두 번째 줄은 "표의 각 행이 정말 내 seed 항체를 가리키는가"를 확인합니다.

## 7. 메시지 읽는 법

모든 명령의 출력은 같은 모양입니다. 줄 하나가 발견 사항 하나이고, 맨 앞
단어가 중요도입니다.

| 표시 | 뜻 | 어떻게 하나 |
|---|---|---|
| **ERROR** | ABD가 거부할 가능성이 큰 문제 | 메시지의 안내대로 고친다. 종료 코드 1 |
| **WARNING** | 사실일 수도 있고 문제일 수도 있는 것 | 읽고 판단한다. 종료 코드는 0 (`--strict`를 붙이면 1) |
| **INFO** | 도구가 한 일의 기록 | 읽지 않아도 된다 |

둘째 단어(`PDB_RESOLUTION`처럼 대문자)는 **finding code**로, 사용자 가이드
[부록 B](user-guide.ko.md#부록-b-자주-만나는-finding-code-빠른-대처표)에
코드별 대처법이 있고, 영어 README에 98개 전부의 설명이 있습니다. 모든
명령은 같은 폴더에 `*.report.json`도 남깁니다. 무엇을 했고 입력·출력
파일의 지문(SHA-256)이 무엇인지 적힌 기록이니, 올린 파일과 함께 보관하세요.

## 8. 자주 막히는 곳

**`python3: command not found`** -- Python이 없거나 이름이 다릅니다. Windows는
`py -3`, 일부 환경은 `python`으로 바꿔 보세요. 2장 참고.

**`No module named abd_prep`** -- 저장소 폴더 밖에서 실행했습니다. `cd`로
`abd-prep` 폴더(안에 `abd_prep/`와 `README.md`가 있는 곳)에 들어가서 다시
실행하세요. 어디서든 쓰고 싶으면 폴더 안에서 `pip install .`을 한 번 하면
`abd-prep quickstart`처럼 짧게 부를 수 있습니다.

**내려받기가 안 된다 (`-> HTTP ...`, `timed out`)** -- 회사 네트워크면 프록시
설정이 필요할 수 있습니다. 터미널에서 `HTTPS_PROXY` 환경 변수를 설정하면
도구가 그대로 따릅니다. 접속하는 곳은 `files.rcsb.org`와 `data.rcsb.org`
둘뿐입니다. 한 번 받은 파일은 `~/.cache/abd-prep/`에 보관되어 다음부터는
네트워크 없이 돌아갑니다.

**`5 file(s) already exist ... pass --force`** -- 같은 폴더에 두 번 실행했습니다.
덮어쓰려면 `--force`, 따로 만들려면 `--out 다른폴더`.

**`WARNING PDB_RESOLUTION`이 신경 쓰인다** -- 1N8Z의 해상도(2.52 Å)는 ABD
권장치(2.5 Å 이하)와 거의 같습니다. 블로그 예제와 같은 구조이니 그대로
써도 되고, 더 선명한 HER2 구조를 쓰고 싶다면
`--preset her2-trastuzumab-6lbx-6b9z`도 있습니다(단, 그 항체는 트라스투주맙의
변이체이며 README에 설명이 있습니다).

**ABD가 파일을 거부했다** -- 거부 메시지를 그대로 들고 사용자 가이드 12장과
13장(FAQ)을 보세요. 세 가지 미확인 규약(체인 ID, FASTA 헤더, CSV 열) 중
하나일 가능성이 가장 큽니다.

## 9. 다음에 읽을 것

- [사용자 가이드](user-guide.ko.md) -- 명령·옵션 전체, PDB/FASTA/CSV 안에
  무엇이 들어 있는지, 메시지별 대처, 용어집.
- [README (영어)](../README.md) -- 기준 명세. finding code 98개 전부.
- [개발자 가이드](developer-guide.ko.md) -- 도구를 고치거나 프리셋을 추가하려면.

## 용어 열 개

| 용어 | 쉬운 뜻 |
|---|---|
| **타겟 항원** | 항체가 붙어야 하는 표적 단백질. 여기서는 HER2 |
| **seed 항체** | 개량·평가의 출발점이 되는 기존 항체. 여기서는 트라스투주맙 |
| **Fv** | 항체에서 항원을 실제로 붙잡는 가변 영역. 중쇄 VH + 경쇄 VL |
| **중쇄 / 경쇄 (H / L)** | 항체를 이루는 두 종류의 단백질 사슬. 긴 것이 중쇄 |
| **PDB** | 단백질 구조(원자 하나하나의 3D 좌표) 파일 형식, 또는 그 공개 데이터베이스 |
| **FASTA** | `>이름` 줄 다음에 아미노산을 글자로 나열한 서열 텍스트 형식 |
| **CSV** | 쉼표로 칸을 나눈 표 텍스트. Excel에서 "CSV UTF-8"로 저장하면 된다 |
| **체인 ID** | PDB 파일 안에서 사슬을 구별하는 한 글자(A, B, H, L ...) |
| **잔기** | 단백질 사슬의 아미노산 한 개. "107잔기" = 아미노산 107개 |
| **프리셋** | "어느 구조에서 어느 체인·범위를 잘라 어떤 이름으로 쓸지"를 미리 정해 둔 묶음 |
