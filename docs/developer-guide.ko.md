# abd-prep 개발자 가이드

`abd-prep`의 내부 구조, 설계 원칙, 각 모듈이 하는 일, 확장·테스트·릴리스
절차를 설명합니다. 도구를 **고치거나 넓히려는** 사람을 위한 문서입니다.
도구를 **쓰는** 방법은 [사용자 가이드](user-guide.ko.md)에, 명령·옵션·finding
code의 기준 명세는 영어 [README](../README.md)에 있습니다.

## 목차

- [1. 개발 환경](#1-개발-환경)
- [2. 설계 원칙](#2-설계-원칙)
- [3. 아키텍처와 데이터 흐름](#3-아키텍처와-데이터-흐름)
- [4. 공통 계층 `abd_prep/common/`](#4-공통-계층-abd_prepcommon)
- [5. PDB 그룹 `abd_prep/pdb/`](#5-pdb-그룹-abd_preppdb)
- [6. FASTA 그룹 `abd_prep/fasta/`](#6-fasta-그룹-abd_prepfasta)
- [7. CSV 그룹 `abd_prep/csvdata/`](#7-csv-그룹-abd_prepcsvdata)
- [8. CLI 계층 `abd_prep/cli.py`](#8-cli-계층-abd_prepclipy)
- [9. report.json 스키마](#9-reportjson-스키마)
- [10. 확장 레시피](#10-확장-레시피)
- [11. 테스트 전략](#11-테스트-전략)
- [12. 릴리스 절차](#12-릴리스-절차)
- [13. 보안 고려 사항](#13-보안-고려-사항)
- [14. 알려진 한계와 후속 과제](#14-알려진-한계와-후속-과제)
- [15. 개발 이력과 주요 결정](#15-개발-이력과-주요-결정)
- [16. 코딩 규약](#16-코딩-규약)

## 1. 개발 환경

필요한 것은 Python 3.8 이상뿐입니다. 서드파티 패키지는 실행에도
테스트에도 없습니다.

```bash
git clone https://github.com/hmkim/abd-prep.git
cd abd-prep
python3 -m abd_prep --help                       # 설치 없이 실행
python3 -m unittest discover -s tests -t .       # 전체 테스트 (네트워크 불필요)
python3 tests/make_seq_fixtures.py --check       # FASTA/CSV 픽스처가 생성기와 바이트 동일한지
```

편집 가능 설치(`pip install -e .`)를 하면 `abd-prep` 명령이 소스 트리를
그대로 가리킵니다. `pytest`가 있다면 `pytest tests/`로도 같은 테스트가
수집됩니다(일부 테스트는 `subTest`를 쓰므로 pytest는 "278 passed, N
subtests"처럼 보고합니다).

**실제 RCSB 파일로 하는 회귀 테스트**(`tests/test_regression_prototype.py`,
`tests/test_split.py`의 프리셋 회귀)는 `~/.cache/abd-prep/`에 원본이 있을 때만
실행되고 없으면 자동으로 건너뜁니다. 한 번 받아 두면 됩니다.

```bash
python3 -m abd_prep pdb fetch 1N8Z 6LBX 6B9Z
```

저장소 구조:

```
abd_prep/
  __init__.py         __version__
  __main__.py         python3 -m abd_prep 진입점
  cli.py              argparse 트리, 액션당 함수 하나, 종료 코드 매핑
  quickstart.py       간단 모드: split + seed FASTA + CSV 예시 + README.txt 를 한 폴더에 (8.1절)
  common/             result(finding 모델), report(report.json), hashing(SHA-256/크기/UTF-8), textfile(선검사), residues(아미노산 표)
  pdb/                pdbfile(고정 칸 레코드, TER/END), ranges(선택 문법), clean, validate, fetch, presets, split
  fasta/              fastafile(파서·정규화·쓰기), roles(역할 토큰·ABD 헤더), validate, crosscheck(F18)
  csvdata/            schema(표준 열·결측 토큰·--map), csvfile(구분자·헤더·행), validate, crosscheck(F22)
tests/
  test_*.py           unittest 스위트
  fixtures/           PDB 3종 + FASTA 7종 + CSV 3종
  make_fixtures.py    PDB 픽스처 생성기 (캐시 원본에서 잘라냄)
  make_seq_fixtures.py FASTA/CSV 픽스처 생성기 (--check 모드)
  helpers.py          임시 디렉터리·CLI 실행 도우미
docs/                 이 문서, 사용자 가이드, 쉬운 사용 방법(튜토리얼)
```

`csvdata/`라는 이름은 표준 라이브러리 `csv`와의 충돌을 피하기 위한
것입니다(패키지 안에서 `import csv`가 자기 자신을 가리키는 사고 방지).
같은 이유로 선택 문법 모듈은 `select.py`가 아니라 `ranges.py`입니다.

## 2. 설계 원칙

코드를 읽거나 바꿀 때 기준이 되는 원칙입니다. 새 기능은 이 원칙과
충돌하지 않아야 합니다.

1. **표준 라이브러리만.** `argparse`, `csv`, `re`, `json`, `hashlib`,
   `urllib.request`, `dataclasses`, `typing`. 사용자가 설치 없이 돌릴 수
   있어야 하고, 보안 검토 대상이 Python 자체로 끝나야 합니다. 서드파티가
   필요한 기능(mmCIF, `.xlsx`, ANARCI)은 선택 extra로만 검토합니다(14장).

2. **형식은 고치고 과학은 바꾸지 않는다.** 잔기 번호·이름, 서열 내용,
   측정값은 절대 바꾸지 않습니다. 바꿀 수 있는 것은 포장(인코딩, 줄바꿈,
   구분자, 헤더, altloc 칸, TER/END, 체인 ID, 변형 잔기의 **명시적** 표준
   변환)뿐입니다. 판단이 필요한 것은 finding으로 보고하고 사용자에게
   맡깁니다. 예: 서열의 `X`를 자동으로 바꾸지 않고 `--replace-ambiguous`를
   요구, 끝 `*`는 지우되 중간 `*`는 오류.

3. **finding code는 고정 식별자, 종료 코드는 계약.** 모든 관찰은
   `Finding(level, code, message, where)`로 기록됩니다. `code`는
   `GROUP_SUBJECT_STATE` 모양의 대문자 상수이고 README 표·테스트·리포트가
   같은 문자열을 공유합니다. 메시지는 바꿔도 코드는 바꾸지 않습니다. 종료
   코드는 0(통과) / 1(검증 오류, `--strict`면 경고도) / 2(입력·옵션 오류)의
   세 가지만 있고, 2는 `InputError` 예외로만 발생합니다.

4. **원문 줄을 칸 단위로만 고친다(PDB).** `Atom`은 80칸 원문 줄을 그대로
   들고 다니며 특정 칸(17 altloc, 22 체인, 7-11 일련번호 ...)만 잘라
   붙입니다. 좌표·B-factor·원자 순서는 건드리지 않으므로 프로토타입 산출물과
   **바이트 단위 동일성**이 구조적으로 보장됩니다. 파싱 후 재조립하는
   방식으로 바꾸면 이 보장이 깨집니다.

5. **추측하지 않는다.** 중쇄/경쇄 역할은 헤더 토큰·옵션으로만 정하고 서열
   내용으로 추정하지 않습니다. framework도 레코드 수로 추정하지 않습니다.
   형식 그룹(`pdb`/`fasta`/`csv`)도 확장자로 자동 판별하지 않고 사용자가
   명령으로 고릅니다. 틀린 추측은 조용히 틀린 파일을 만들기 때문입니다.

6. **미확인 외부 규약은 상수 한 곳에 격리한다.** ABD 콘솔에서 글자 단위로
   확인되지 않은 규약 3가지(PDB 체인 라벨, FASTA 헤더, CSV 열 스키마)는 각각
   `pdb/presets.py`(`Part.out_chain`), `fasta/roles.py`(`ABD_HEADER_ID`,
   `fasta_header`), `csvdata/schema.py`에만 있습니다. 확인되면 그 파일과
   README만 고칩니다.

7. **테스트는 네트워크 없이, 픽스처는 출처가 명시된 공개 데이터만.**
   네트워크 호출은 `opener` 주입으로 테스트에서 가짜로 바꿉니다. 픽스처는
   RCSB 공개 항목에서 잘라 낸 것(REMARK/헤더에 출처 표기), 명시된 합성
   서열, 지어낸 숫자만 담습니다. 실측 데이터는 넣지 않습니다.

## 3. 아키텍처와 데이터 흐름

세 형식 그룹은 서로 독립이고 `common/`만 공유합니다. 한 명령의 흐름은
어느 그룹이든 같은 모양입니다.

```
argv
 │  cli.build_parser()          argparse: GROUP -> ACTION -> 옵션
 ▼
cmd_<group>_<action>(args, argv) -> int
 │
 ├─ 입력 존재 확인 (_require_file)            → InputError → 종료 2
 ├─ Report(command=render_command(argv))      리포트 골격 (명령·tool 버전·시각)
 ├─ 선검사  textfile.precheck()               크기 100MB / 인코딩 / BOM / 확장자   → FILE_* findings
 ├─ 파싱    pdbfile.read_pdb | fastafile.read_fasta | csvfile.read_layout+iter_rows
 ├─ 처리    clean_atoms | validate_atoms | validate_records | crosscheck | process_csv
 │            └─ ValidationResult 에 finding 축적, details dict 에 수치
 ├─ 쓰기    write_pdb | write_fasta | process_csv(out_path=...)   (clean/split 만)
 ├─ report.add_input/add_output  (경로·SHA-256·바이트)
 └─ _emit()  텍스트 표 또는 --json, report.json 저장
 ▼
result.exit_code(strict)  → 0 / 1
```

`ValidationResult`는 그룹 사이를 오가는 유일한 "결과" 타입입니다. 라이브러리
함수들은 `result`를 인자로 받아 finding을 **덧붙이고**, 수치는 별도 `details`
dict로 돌려줍니다. 출력(표·JSON·파일)은 `cli.py`에서만 일어나고, 라이브러리
모듈은 `print`하지 않습니다.

라이브러리로 쓸 때의 최소 예:

```python
from abd_prep.common.result import ValidationResult
from abd_prep.pdb.clean import CleanOptions, clean_file
from abd_prep.pdb.ranges import ChainSelection, parse_rename
from abd_prep.pdb.validate import validate_file

opts = CleanOptions(selection=ChainSelection.parse(["A:1-107", "B:1-120"]),
                    rename=parse_rename("A=L,B=H"))
outcome = clean_file("1N8Z.pdb", "fv.pdb", opts)      # CleanOutcome(atoms, result, details, ters)
print(outcome.result.format_table("[clean]"))

res, details = validate_file("fv.pdb")
print(res.exit_code(strict=False), details["chains"]["H"]["residues"])
```

## 4. 공통 계층 `abd_prep/common/`

### 4.1 `result.py` -- finding 모델과 종료 코드

```python
LEVELS = ("info", "warning", "error")
EXIT_OK, EXIT_VALIDATION, EXIT_INPUT = 0, 1, 2

class InputError(Exception): ...          # CLI 가 종료 2 로 매핑

@dataclass
class Finding:
    level: str; code: str; message: str; where: str = ""

@dataclass
class ValidationResult:
    findings: List[Finding]
    def info(code, message, where="") / warning(...) / error(...)
    def extend(other)                     # 다른 결과 합치기
    def errors / warnings / has_errors / has_warnings / codes
    def exit_code(strict=False) -> int    # error→1, strict and warning→1, else 0
    def to_list() / to_json()
    def format_table(title) -> str        # CLI 텍스트 표
```

`where`는 사람이 읽는 위치 문자열입니다(`file`, `A`, `C:580-591`,
`H:B<-H`, `line 8486`, `header`, `rows`, `value` ...). 파일이 여럿일 때 CLI가
앞에 파일 이름을 붙입니다. `format_table`은 코드·위치 열 너비를 내용에 맞춰
정렬하고 끝에 `summary: N error(s), N warning(s), N info`를 씁니다.

### 4.2 `report.py` -- report.json

```python
Report(command: str, path_mode: "basename"|"full")
  .add_input(kind, path)   # kind: pdb | fasta | csv ; sha256 + bytes 를 즉시 계산
  .add_output(kind, path)
  .result: ValidationResult
  .details: dict
  .to_dict() / .to_json() / .write(path)
render_command(argv, path_mode)   # "abd-prep pdb clean in.pdb out.pdb ..." (경로는 mode 에 따라 축약)
```

`generated_at`은 UTC ISO 8601, `tool`은 이름·`__version__`·Python 버전입니다.
경로 기본값이 basename인 이유는 리포트가 메일·티켓에 첨부될 때 홈 디렉터리와
사용자 이름이 새지 않게 하려는 것입니다.

### 4.3 `hashing.py` -- 크기·SHA-256·UTF-8

- `sha256_file(path)`: 1 MiB 청크 스트리밍.
- `ABD_MAX_BYTES = 100 MiB`, `exceeds_abd_limit(path)`.
- `check_utf8(path, encoding="utf-8") -> (ok, message)`: 청크 단위로 디코딩하되
  청크 경계에서 잘린 다중 바이트 문자를 틀리게 보고하지 않도록
  `codecs.getincrementaldecoder`를 씁니다. 오류 위치는 `offset + exc.start`
  (파일 처음부터의 바이트 위치).
- `has_utf8_bom(path)`.

### 4.4 `textfile.py` -- FASTA/CSV 공통 선검사

`precheck(path, encoding, result, action) -> Optional[str]`가 한 번에
처리합니다. 돌려주는 값은 실제로 쓸 코덱 이름이고, `None`이면 파일을 읽을 수
없다는 뜻(오류 finding은 이미 기록됨)입니다.

| 단계 | finding |
|---|---|
| 알 수 없는 코덱 이름 (`--encoding foo`) | `InputError` (종료 2) |
| 100 MB 이상 | `FILE_TOO_LARGE` 오류 (계속 진행) |
| 디코딩 실패 | UTF-8이면 `FILE_NOT_UTF8`, 지정 코덱이면 `FILE_ENCODING` -- 오류, `None` 반환 |
| 지정 코덱으로 읽혔는데 ASCII뿐 | `FILE_ENCODING_ASCII` 정보 (이미 UTF-8) |
| 지정 코덱으로 읽힘 (validate) | `FILE_ENCODING_INPUT` 경고 -- ABD는 UTF-8 필요 |
| 지정 코덱으로 읽힘 (clean) | `FILE_ENCODING_CONVERTED` 정보 -- 출력은 UTF-8 |
| UTF-8 BOM | `FILE_BOM` 정보; 반환 코덱이 `utf-8-sig`로 바뀌어 BOM이 데이터로 들어오지 않음 |

`check_extension(path, allowed, group, result)`: `.xlsx`/`.xls`/`.xlsm`/`.xlsb`는
"Excel에서 CSV UTF-8로 저장" 안내와 함께 `InputError`, 그 외 낯선 확장자는
`FILE_EXTENSION` 경고만 하고 파싱은 계속합니다(`.faa`, `.txt` 등이 흔하기
때문). PDB 그룹은 자체 `PDB_EXTENSION` 경고를 씁니다.

### 4.5 `residues.py` -- 아미노산 표

- `AA3_TO_1` / `AA1_TO_3`, `STANDARD_AA3`, `STANDARD_AA1`, `AMBIGUOUS_AA1 =
  {X,B,Z,J,U,O}`.
- `ModifiedResidue(NamedTuple)` 변환표: `MSE→MET`(SE→SD), `CSO/CSD/OCS/CME→CYS`,
  `SEP→SER`, `TPO→THR`, `PTR→TYR`, `HYP→PRO`, `MLY/KCX→LYS`, `PCA→GLN`(OE→OE1).
  변환은 "원자 이름 바꾸기 → 표준 잔기의 정규 원자 집합(`standard_atoms`)에
  없는 원자 제거"라는 **일반 규칙**으로 수행되므로 잔기별 삭제 목록을 따로
  관리하지 않습니다. 새 변형 잔기를 추가할 때도 모체와 이름 변경 쌍만
  적으면 됩니다(10장).
- `check_protein_alphabet(seq) -> (bad_positions, looks_nucleic)`: FASTA와 CSV
  서열 열이 공유.
- `element_for_atom`, `is_hydrogen_name`: 원소 칸이 비어 있을 때의 폴백.
- `sequence_from_residues(resnames)`: 3글자 → 1글자, 모르면 `X`.

## 5. PDB 그룹 `abd_prep/pdb/`

### 5.1 `pdbfile.py` -- 고정 칸 레코드

`Atom`은 원문 한 줄(`line`)과 `model`(0 = MODEL 블록 밖), `lineno`를 갖는
데이터클래스이고, 모든 필드는 **칸을 잘라 읽는 프로퍼티**입니다
(0-기준 슬라이스; PDB 명세의 1-기준 칸 번호와 비교해 두세요).

| 필드 | 슬라이스 | PDB 칸 | 비고 |
|---|---|---|---|
| `record` | `[0:6]` | 1-6 | `ATOM` / `HETATM` |
| `serial` | `[6:11]` | 7-11 | 정수가 아니면 `None` |
| `name` / `name_field` | `[12:16]` | 13-16 | `name_field`는 공백 포함 원문(altloc 키에 사용) |
| `altloc` | `[16]` | 17 | 공백이면 대체 위치 없음 |
| `resname` | `[17:20]` | 18-20 | |
| `chain` | `[21]` | 22 | |
| `resseq` | `[22:26]` | 23-26 | 정수 아니면 `InputError`(hybrid-36은 범위 밖) |
| `icode` | `[26]` | 27 | 삽입 코드, 없으면 `" "` |
| `xyz` | `[30:38]`,`[38:46]`,`[46:54]` | 31-54 | |
| `occupancy` | `[54:60]` | 55-60 | 없으면 `None` |
| `element` | `[76:78]` | 77-78 | |

편집 메서드(`set_altloc`, `set_chain`, `set_record`, `set_resname`, `set_name`,
`set_element`, `set_serial`)는 해당 칸만 `_splice`로 바꿉니다. 짧은 줄은 80칸으로
패딩한 뒤 바꿉니다. `is_hydrogen`은 원소 칸 → 이름 폴백 순서로 판정하며
`clean --strip-h`와 `validate`가 **같은 프로퍼티**를 써서 두 명령의 수소
집계가 항상 일치합니다.

파일 입출력은 `FILE_ENCODING = "latin-1"`로 합니다. PDB는 ASCII이지만
간혹 섞여 들어오는 비ASCII 바이트(편집기 유니코드 문자 등)를 latin-1은
1바이트 = 1문자로 **손실 없이 왕복**시키므로 트레이스백 없이 그대로 쓰고
`PDB_NON_ASCII` 경고만 냅니다.

`read_pdb(path) -> PdbFile`은 `ATOM`/`HETATM`을 `Atom`으로, `TER`을
`TerRecord`로 모으고 `MODEL`/`ENDMDL`로 모델 번호를 매기며, `END` 유무를
기록합니다. 다른 레코드는 버립니다. `iter_residues(atoms)`는 파일 순서로
`(residue_key, [atoms])`를 냅니다.

`write_pdb(atoms, out, renumber) -> ters`는 (체인 ID 치환은 `clean_atoms`가
이미 끝낸) 원자 줄을 쓰면서 선택적으로 일련번호를 1부터 다시 매기고,
**체인이 바뀌는 지점마다 `TER`을 새로 만들어** 직전 잔기(이름·체인·번호·삽입
코드)를 적습니다. 원본 `TER`은 쓰지 않습니다. 마지막에 `END`.
`ter_serial_collisions(atoms, ters)`는 재생성한 TER 일련번호가 원자
일련번호와 겹치는지(원본에 체인 사이 TER이 없던 경우) 알려 주어 `clean_file`이
`PDB_SERIAL_DUPLICATE` 경고 + `--renumber` 안내를 냅니다.

### 5.2 `ranges.py` -- 선택 문법

```
SPEC := CHAIN | CHAIN ":" RANGE ("," RANGE)*
RANGE := RES | RES "-" RES
RES := [-]digits [icode-letter]
```

- `ResidueRange(lo, lo_icode, hi, hi_icode, hi_exact)`: 상한에 삽입 코드를
  **쓰지 않으면** `hi_exact=False`로 그 번호의 모든 삽입 잔기를 포함하고
  (`1-82` ⊇ 82, 82A, 82B...), **쓰면** 정확히 그 잔기까지입니다. 삽입 코드
  순서는 공백 < `A` < `B` < ... 이며 내부적으로 `_ICODE_MAX = "\x7f"`를
  "모든 코드 포함" 상한으로 씁니다.
- `ChainSelection.parse(specs)`: 체인별 범위 목록. `wants_chain(chain)`,
  `wants(chain, resseq, icode)`. 비어 있으면 모든 체인 선택. `describe()`는
  리포트용 `{chain: "spec"}`.
- `parse_rename("A=L,B=H") -> {"A": "L", "B": "H"}`. ID는 한 글자.
- 문법 오류는 `RangeSyntaxError`로, CLI가 `InputError`(종료 2)로 바꿉니다.
  `10-5`는 "start is after end".

### 5.3 `clean.py` -- 정리 알고리즘

`clean_atoms(pdb, opts) -> CleanOutcome(atoms, result, details, ters)`의 순서:

1. **첫 모델만**: `atom.model in (0, 첫 MODEL 번호)`. 둘 이상이면
   `PDB_MULTI_MODEL` 경고.
2. **잔기 수준 이질성 해소** -- `residue_winners(atoms, candidate)`: 같은
   `(chain, resseq, icode)`에 잔기 이름이 둘 이상이면(6B9Z B175 `ACSO`/`BCYS`)
   컨포머별 점유율 합이 큰 **이름**을 먼저 고릅니다. 이때 `candidate`
   함수(`survives_het_filter`)로 HETATM 필터를 통과하지 못할 컨포머는 경쟁에서
   제외해, 지워질 CSO가 살아남을 CYS를 밀어내는 일을 막습니다. 결과는
   `PDB_RESIDUE_HETEROGENEITY` 정보, 제거된 변형 컨포머는
   `PDB_MODIFIED_CONFORMER_REMOVED` 정보(구멍이 아니므로 경고가 아님).
3. **원자 수준 altloc 선택** -- `select_altlocs(atoms, candidate)`: 키
   `(chain, resseq text, icode, name_field)`마다 점유율 최대 컨포머(동률이면
   파일에서 먼저 나온 것). 이 알고리즘과 동률 규칙은 프로토타입과 같아 참조
   산출물 4개가 바이트 동일합니다. 선택된 원자는 `set_altloc(" ")`.
   `PDB_ALTLOC_RESOLVED` 정보에 영향 원자·잔기·버린 원자 수.
4. **선택 적용** -- `ChainSelection.wants`. 지정한 체인이 파일에 없으면
   `PDB_CHAIN_NOT_FOUND`, 결과가 비면 `PDB_SELECTION_EMPTY`(HETATM만 걸려
   지워진 경우는 `--keep-het` 힌트), 원자가 하나도 없으면 `PDB_NO_ATOMS`.
5. **HETATM 처리** -- `polymer_like_hetero_residues(atoms)`가 "백본 삼총사
   `N`/`CA`/`C`를 모두 갖고, 같은 체인에 `ATOM` 레코드가 있는 HETATM 잔기"를
   폴리머형으로 분류합니다(물·글리칸·이온은 탈락, MSE/CSO/PCA는 통과).
   - 기본: 모든 HETATM 제거. 폴리머형이 선택 범위 안에서 지워지면
     `PDB_MODIFIED_RESIDUE_REMOVED` 경고(구멍), 변환 규칙이 없는 것은
     `PDB_MODIFIED_RESIDUE_NO_RULE` 경고.
   - `--std-residues`: 변환 규칙이 있는 폴리머형은 `residues.ModifiedResidue`에
     따라 `ATOM`으로 바꾸고 이름 변경 후 정규 원자 집합 밖 원자를 제거
     (`PDB_MODIFIED_RESIDUE_CONVERTED` 정보에 버린 원자 목록).
   - `--keep-het`: 전부 유지.
6. **수소 제거**(`--strip-h`): `Atom.is_hydrogen`. `PDB_HYDROGENS_REMOVED` 정보.
7. **체인 이름 변경**(`--rename`): 마지막 단계에서 원자 줄마다 `set_chain`으로
   치환합니다. 각 줄의 **원본** 체인 ID를 보고 바꾸므로 `A=B,B=A` 교환도
   안전합니다. 두 원본 체인이 한 출력 ID로 합쳐지면 `PDB_RENAME_COLLISION`
   오류(출력은 그래도 씀, 종료 1).

`clean_file(in, out, opts)`는 `read_pdb → clean_atoms → write_pdb`를 묶고
TER 일련번호 충돌 경고를 추가합니다. `details`에는 `atoms_in/out`,
`hetatm_removed{resname: n}`, `rename`, `chains_out{chain: residues/atoms/first/last}`,
`ter_records_written`, 변환·제거 목록이 들어갑니다.

### 5.4 `validate.py` -- 검증 규칙

`validate_atoms(pdb, opts) -> (ValidationResult, details)`:

| 검사 | code | 수준 |
|---|---|---|
| 원자 없음 | `PDB_NO_ATOMS` | 오류 |
| 17번 칸이 빈 칸이 아닌 원자 | `PDB_ALTLOC_PRESENT` | 오류 |
| `END` 없음 | `PDB_END_MISSING` | 오류 |
| 체인 수 ≠ `expect_chains` | `PDB_CHAIN_COUNT` | 오류 |
| `TER` 규칙(`_check_ters`): 경계에 없음 / 경계가 아닌 곳 / 앞 잔기와 다름 / 한 경계에 여럿 | `PDB_TER_MISSING` / `PDB_TER_MISPLACED` / `PDB_TER_MISMATCH` / `PDB_TER_EXTRA` | 오류 |
| HETATM 존재 | `PDB_HETATM_PRESENT` | 경고 (`allow_het`면 정보) |
| 연속 잔기 CA-CA > `CA_GAP_ANGSTROM = 4.2` | `PDB_GAP_CA_DISTANCE` | 경고 (`expected_gaps`에 있으면 정보 + "expected") |
| 잔기 번호 점프·감소 | `PDB_GAP_NUMBERING` | 정보 |
| `ATOM`에 비표준 잔기 이름 | `PDB_NONSTANDARD_RESNAME` | 경고 |
| 체인 H/L 잔기 수가 `FV_RESIDUE_RANGE` (L 105-115, H 110-125) 밖, 또는 `expect_residues` 불일치 | `PDB_RESIDUE_COUNT` | 경고 |
| 같은 번호에 다른 잔기 이름 | `PDB_RESIDUE_DUPLICATE` | 경고 |
| 모델 여럿 / 빈 체인 ID / 일련번호 중복 / 비ASCII | `PDB_MULTI_MODEL` / `PDB_CHAIN_BLANK` / `PDB_SERIAL_DUPLICATE` / `PDB_NON_ASCII` | 경고 |
| 수소 존재 | `PDB_HYDROGENS_PRESENT` | 정보 |

`ValidateOptions(allow_het, expect_chains, expect_residues{chain: (lo, hi)},
expected_gaps{(chain, before, after)})`는 프리셋이 채워 줍니다.
`summarize_chains`가 체인별 `ChainSummary`(first/last/residues/atoms/sequence/
nonstandard)를 만들고 `inspect_atoms`/`format_inspect`/`to_fasta`가 그것을
표·FASTA로 냅니다. `to_fasta`는 `fasta.roles.fasta_header`를 써서 체인 `H`/`L`을
ABD 헤더로 적습니다 -- 헤더 규약이 바뀌어도 한 곳만 고치면 되는 이유입니다.

### 5.5 `fetch.py` -- RCSB 다운로드

- 상수: `RCSB_FILES_URL = https://files.rcsb.org/download/`,
  `RCSB_ENTRY_URL = https://data.rcsb.org/rest/v1/core/entry/`,
  `DEFAULT_CACHE_DIR = ~/.cache/abd-prep`, `CACHE_ENV = ABD_PREP_CACHE`,
  `TIMEOUT_SECONDS = 60`, `MAX_DOWNLOAD_BYTES = 512 MiB`.
- `normalize_pdb_id(text)`: `^[1-9][A-Za-z0-9]{3}$` 검증 후 대문자화. 이 검증을
  통과한 문자열만 URL과 캐시 **경로**에 쓰이므로 `../x` 같은 경로 이탈이
  불가능합니다.
- `fetch_pdb(pdb_id, cache, refresh, opener)`: 캐시 `<ID>.pdb`가 있으면 그것,
  없으면 `.part`로 받아 `os.replace`(부분 파일이 캐시에 남지 않음). 응답 머리가
  PDB 레코드로 시작하는지 확인합니다. `.pdb`가 404인데 entry API는 200이면
  "mmCIF 전용 항목", 둘 다 404면 "없는 항목"으로 `FetchError`(=`InputError`)
  메시지를 나눕니다.
- `fetch_entry(pdb_id, ...)`: entry JSON에서 제목·개정일·공개일·해상도·실험
  방법·모델 수만 추려 `EntryMeta`로, `<ID>.entry.json`에 캐시. 실패해도 파일
  처리는 계속되고(`meta_error`), `split`이 `PDB_ENTRY_META_UNAVAILABLE` 정보를
  냅니다.
- `opener: Callable[[url], (status, bytes)]`를 주입할 수 있어 테스트는
  네트워크 없이 가짜 응답으로 모든 분기를 돕니다(`tests/test_fetch.py`).

### 5.6 `presets.py` -- 프리셋

```python
@dataclass(frozen=True)
class Part:
    role: str                 # "target" | "heavy" | "light"
    pdb_id: str               # RCSB ID, 또는 로컬 파일이면 "input"
    selection: str            # "C:489-607" -- 체인 하나, 범위 여러 개 가능
    out_chain: Optional[str]  # None 이면 원본 체인 ID 유지

@dataclass(frozen=True)
class Preset:
    name, description, parts: Tuple[Part, ...]
    expect_residues: Dict[str, int]                 # 출력 체인 -> 측정된 잔기 수
    expected_gaps: FrozenSet[Tuple[str, int, int]]  # (출력 체인, 앞 잔기, 뒤 잔기)
    notes: Tuple[str, ...]
```

`PRESETS` dict에 `her2-trastuzumab`(기본)과 `her2-trastuzumab-6lbx-6b9z`가
있습니다. `validate_parts`는 파트가 하나도 없는 경우, 역할 중복, heavy/light가
서로 다른 항목에서 오는 경우(한 파일에 쓰이므로 같은 좌표계여야 함),
heavy/light가 같은 원본 체인을 고르는 경우, 두 파트가 같은 출력 체인 ID로
쓰이는 경우를 거부하고, 모듈 import 시점에 `PRESETS` 전체에 대해 한 번
실행됩니다. `output_groups(parts)`는 `{"target": [...], "fv": [...]}`로 묶고
`output_filename(name, group)`이 `<name>_target.pdb` / `<name>_fv.pdb`를 만듭니다.
`cli._preset_options_for`는 검사 대상 파일의 체인 ID 집합을 각 그룹의 출력
체인 집합과 비교해 어느 파트인지 알아내고 `ValidateOptions`를 채웁니다.

### 5.7 `split.py` -- fetch + clean + validate

`run_split(opts, opener, command) -> SplitOutcome(written, report, report_path)`:

1. `resolve_sources`: 파트들이 가리키는 항목을 `fetch`(또는 `--input` 로컬
   파일)로 모아 `details["sources"]`에 항목 메타(개정일·해상도·방법)를 기록.
   해상도가 2.5 Å보다 거칠면 `PDB_RESOLUTION` 경고.
2. 그룹(`target`, `fv`)마다 해당 파트들의 선택을 합쳐 `clean_atoms` →
   `write_pdb` → `validate_file`(프리셋 기대치 적용). finding `where`에
   `target/` 또는 `fv/` 접두를 붙입니다(`_prefix`).
3. `report.json`을 `<name>.report.json`에 저장(`report_path=""`면 생략).

`derive_name(parts, input_path)`은 `--name`이 없을 때 프리셋 이름, 항목 ID,
또는 입력 파일 이름(확장자 제외)에서 출력 이름을 만들고, `run_split`은 이름
(`--name` 포함)에 경로 구분자가 있거나 `.`/`..`이면 `InputError`로 거부합니다.

## 6. FASTA 그룹 `abd_prep/fasta/`

### 6.1 `fastafile.py` -- 파서·정규화·쓰기

```python
@dataclass
class Record:
    id: str; description: str; sequence: str; lineno: int
    role: Optional[str] = None          # "H" | "L" | "HL" | None
    role_source: Optional[str] = None   # "option" | "framework" | "header"
    header -> str ; length -> int ; to_dict()

@dataclass
class FastaFile:
    records: List[Record]; path: str
    by_id() -> Dict[str, Record]
```

`read_fasta(source, result, path="", encoding=...)`는 파일 경로·텍스트
스트림·줄 목록을 모두 받습니다. 헤더 줄은 `>` 다음을 `split(None, 1)`로
나눠 첫 토큰이 ID, 나머지가 설명입니다(탭도 구분자). 서열 정규화(F13)는
**형식만** 건드립니다.

| 처리 | finding | 비고 |
|---|---|---|
| CR/LF 감지 | `FASTA_CRLF` 정보 | 쓰기는 항상 LF |
| 공백·숫자 제거 | `FASTA_WHITESPACE_REMOVED` 정보 | 번호 붙은 정렬 출력 대응 |
| `-`·`.` 제거 | `FASTA_GAPS_REMOVED` 정보 | 정렬 갭 |
| 소문자 → 대문자 | `FASTA_LOWERCASE` 정보 | |
| 끝의 `*` 전부 제거 | `FASTA_TRAILING_STOP` 정보 | 개수 보고 |
| 중간의 `*` | `FASTA_INTERNAL_STOP` **오류** | 제거하지 않음; 위치는 실제 서열 기준 |
| 첫 `>` 앞의 서열 텍스트 | `FASTA_NO_HEADER` 오류 | |
| `>`가 하나도 없음 / ID 없는 헤더 / 서열 없는 헤더 | `FASTA_NO_RECORDS` / `FASTA_EMPTY_HEADER` / `FASTA_EMPTY_SEQUENCE` 오류 | |
| `MAX_SEQUENCE_LENGTH = 10_000` 초과 | `FASTA_TOO_LONG` 오류 | 항체 구조물이 아님 |
| 같은 ID 반복 | `FASTA_DUPLICATE_ID` 경고 | |

`format_record(record, wrap, rewrite_header)`와 `write_fasta(records, out,
wrap)`는 UTF-8·LF로 쓰고, `rewrite_header=True`(기본)면 `roles.fasta_header`로
`>H <원래 헤더>` 모양을 만듭니다. 서열 글자는 정규화 뒤 그대로입니다.

### 6.2 `roles.py` -- 역할 추정과 ABD 헤더

이 파일에는 다른 곳에 없는 두 가지가 있습니다.

1. **ABD 헤더 규약** `ABD_HEADER_ID = {"H": "H", "L": "L", "HL": "scFv"}`와
   `fasta_header(role_or_chain, description)`. `pdb inspect --fasta`,
   `fasta clean`, 테스트가 모두 이것을 읽습니다. 콘솔에서 규약이 확인되면
   **이 dict(필요하면 `fasta_header`)만** 고칩니다.
2. **토큰 사전** `HEAVY_TOKENS`, `LIGHT_TOKENS`, `HEAVY_PREFIXES`(`ighv`, `ighg`,
   ...), `LIGHT_PREFIXES`(`igkv`, `iglv`, ...), `SCFV_TOKENS`. `header_tokens`는
   헤더를 `[^A-Za-z0-9]+`로 쪼개 소문자화하므로 RCSB의
   `1N8Z_2|Chains B|... heavy chain` 같은 헤더도 처리됩니다.

`infer_role(header) -> (role, matched_tokens)`: scFv 토큰이 있으면 `HL`,
heavy·light 토큰이 함께 있으면 `None`(모호), 한쪽만 있으면 그 역할, 없으면
`None`. `assign_roles(records, result, framework, heavy_id, light_id)`가 우선순위
(옵션 → 단일 레코드 framework → 헤더)로 `Record.role`/`role_source`를 채우고
`FASTA_ROLE_UNKNOWN`(igg 오류·그 외 경고), `FASTA_ROLE_AMBIGUOUS`,
`FASTA_ROLE_CONFLICT`(nanobody에 light 토큰), `FASTA_ROLE_OVERRIDDEN`(헤더와 옵션이
다름, 정보)을 냅니다. `check_framework`는 레코드 수·역할 조합을
`FRAMEWORKS = ("igg", "nanobody", "scfv")` 규칙과 비교해
`FASTA_FRAMEWORK_MISMATCH`를, framework가 없으면 `FASTA_FRAMEWORK_UNSPECIFIED`
정보를 냅니다. 서열 내용으로 역할을 추정하는 코드는 **의도적으로 없습니다**.

### 6.3 `validate.py` -- 알파벳·정합·휴리스틱

- `check_alphabet(rec, result, replace)`: `residues.check_protein_alphabet`로
  비표준 글자 위치를 모아 `FASTA_BAD_ALPHABET` 오류(`X@10, B@50` 식, 최대 8개
  표시). `--replace-ambiguous` 매핑이 있으면 그 글자를 치환하고
  `FASTA_AMBIGUOUS_REPLACED` 정보 -- `parse_replacements`는 영문자→영문자만
  허용하고 `*`·기호·숫자는 `InputError`. A/C/G/T/U/N만으로 된
  `NUCLEIC_MIN_LENGTH = 20` 이상 서열은 `FASTA_NUCLEIC` 오류(짧은 펩타이드
  오판 방지).
- `check_heuristics(records, framework, result)` (`--no-heuristics`로 끔,
  `FASTA_HEURISTICS_OFF` 정보): `LENGTH_RANGE` = VL 95-130, VH 105-140,
  VHH 100-140, scFv 220-280 → `FASTA_LENGTH_RANGE` 경고. 구조 쪽
  `FV_RESIDUE_RANGE`(105-115/110-125)보다 넓은 이유는 seed 서열이 몇 잔기의
  여분을 정당하게 가질 수 있기 때문입니다. `looks_like_signal_peptide`: `M`으로
  시작하고 앞 `SIGNAL_SCAN = 30`잔기 안에 길이 `SIGNAL_WINDOW = 8` 창의 소수성
  (`AILMFVWC`) 글자가 `SIGNAL_MIN_HYDROPHOBIC = 7`개 이상인 창이 있으면
  `FASTA_SIGNAL_PEPTIDE` 경고. 같은 서열 둘은 `FASTA_DUPLICATE_SEQUENCE` 경고.
- `check_scfv_linker`: 역할이 `HL`인 레코드에서 `LINKER_RE = (?:G{3,4}S){2,}`를
  찾아 `FASTA_SCFV_LINKER` 정보(없으면 "no common pattern", 서열은 그대로).
- `validate_records(records, opts, result)`가 위를 순서대로 묶고,
  `inspect_records`/`format_inspect`가 `inspect` 출력을 만듭니다.

### 6.4 `crosscheck.py` -- F18: 같은 분자인가

입력은 FASTA 레코드들과 PDB 파일, 출력은 체인별 `Alignment`입니다.

1. `pdb_chain_residues(path)`: PDB를 읽어 체인마다 `ChainResidues(labels,
   letters, index, gaps)`를 만듭니다. `index`는 **gap 인식 가상 인덱스** --
   파일 순서로 잔기를 걸으며 번호가 1보다 크게 뛰면 그만큼 인덱스를 건너뛰고
   (결손), 삽입 코드(82 → 82A → 82B)는 연속으로 셉니다. 변형 잔기 HETATM은
   모체 아미노산 글자로 포함합니다.
2. `resolve_chain_map(chains, chain_map, pdb_path)`: PDB 체인 ID가 `H`/`L`이면
   그대로 역할에 대응, 아니면 `--chain-map H=B,L=A`가 필수(없으면 `InputError`
   -- 사용자 가이드 8.6절의 메시지).
3. `align(chain, seq)`: 오프셋 `o`를 `-span+1 .. len(seq)-1` 전 범위로
   돌리며 `(matches, -mismatches, -|o|)` 튜플이 최대인 오프셋을 고릅니다
   (겹침이 `min(5, len)` 미만인 오프셋은 제외). 전수 탐색 O(서열 길이 × 잔기
   수)이지만 항체 입력은 ms, 허용 최대(10,000 × 1,000)에서 약 2초라 최적화하지
   않았습니다. 결과 `Alignment(offset, overlap, matches, mismatches[],
   uncovered[], gaps_in_range[], extra_n, extra_c)`.
4. `crosscheck(records, pdb_path, chain_map, result)`가 역할별로 레코드를
   골라(scFv `HL` 레코드는 H·L 두 체인 모두에 대조) finding으로 바꿉니다:

| 조건 | code | 수준 |
|---|---|---|
| 동일성 요약 (항상) | `FASTA_PDB_IDENTITY` | 정보 |
| 구조에 없는 잔기(결손) | `FASTA_PDB_GAP` | 정보 |
| 서열이 구조 밖으로 더 긺 (N/C 말단) | `FASTA_EXTRA_RESIDUES` | 정보 |
| 구조 잔기가 서열 밖 | `FASTA_PDB_UNCOVERED` | 경고 |
| 둘 다 있는데 다른 잔기 (최대 `LIST_LIMIT = 20`개 나열) | `FASTA_PDB_MISMATCH` | 오류 |
| 동일성 < `MIN_IDENTITY = 0.5` 또는 일치 < `MIN_MATCHES = 10` | `FASTA_PDB_NO_CORRESPONDENCE` (불일치 목록 생략) | 오류 |
| 매핑된 역할에 레코드가 없음 | `FASTA_PDB_ROLE_MISSING` | 경고 |
| 비교된 레코드가 하나도 없음 | `FASTA_PDB_NOT_COMPARED` | 오류 |

마지막 항목이 오류인 이유: 사용자가 요청한 검사가 아무것도 비교하지 않고
종료 0으로 끝나면 통과한 것처럼 보이기 때문입니다.

## 7. CSV 그룹 `abd_prep/csvdata/`

### 7.1 `schema.py` -- 표준 열과 `--map`

**ABD의 실제 CSV 스키마가 확인되면 이 파일만 고칩니다.**

- `STANDARD_COLUMNS = (id, sequence, heavy, light, value, unit, target, assay,
  notes)` -- 출력 순서이기도 함. `REQUIRED_COLUMNS = (id,)`.
  `SEQUENCE_GROUPS = ((sequence,), (heavy, light))` 중 하나가 있으면 "서열 있음".
  `NUMERIC_COLUMNS = (value,)`.
- `MISSING_TOKENS = {"", na, n/a, nan, none, null, -, --, nd, n.d.}` (공백
  제거·소문자 비교). README 표와 validate 규칙이 모두 이 상수를 참조합니다.
- `classify_number(value) -> (kind, float|None)`: `number` / `missing` /
  `decimal_comma`(`^\d+,\d+$`) / `censored`(`<0.1`, `>=100`) / `text`.
- `parse_map(["id=Clone,sequence=AA_seq", "value=KD, nM"]) -> {std: header}`:
  반복 가능, 한 `--map` 안에서는 쉼표로 분리하므로 헤더에 쉼표가 있으면
  단독 `--map`으로 줍니다. 모르는 표준 열 이름은 `InputError`.
- `resolve_columns(header, mapping) -> (std->index, unmapped indexes)`: 매핑에
  없는 헤더가 표준 이름과 (대소문자 무시) 같으면 자동 매핑. 매핑이 가리킨
  헤더가 파일에 없으면 실제 헤더 목록을 담은 `InputError`.

### 7.2 `csvfile.py` -- 구분자·헤더·행

- `parse_delimiter("tab" | "\\t" | ";" ...)`: `--delimiter` 옵션 해석.
- `sniff_delimiter(sample)`: 앞 `SNIFF_BYTES = 64 KiB`를 `csv.Sniffer`에
  `SNIFF_DELIMITERS = ",;\t|"`로 넘김. 실패하면 `read_layout`이 헤더 줄에 후보
  구분자가 **하나만** 들어 있는지 보고 그것을 채택, 그래도 모르면 `,` +
  `CSV_DELIMITER_GUESSED` 경고. 결과는 `CSV_DELIMITER` 정보로 기록.
- `read_layout(path, codec, delimiter, result) -> CsvLayout(delimiter, header,
  header_lineno)`: 헤더 없음 / 빈 이름 / 중복 이름 검사
  (`CSV_NO_HEADER` / `CSV_HEADER_EMPTY_NAME` / `CSV_HEADER_DUPLICATE` 오류).
- `iter_rows(path, codec, layout, stats)`: `csv.reader`로 스트리밍하며
  `(lineno, cells)`를 냅니다. 행 번호는 `reader.line_num`이라 따옴표 안의
  줄바꿈이 있어도 정확합니다. 빈 행은 건너뛰고 `RowStats`에 셉니다
  (`CSV_BLANK_ROWS` 정보). `csv.field_size_limit`은 ABD 파일 상한(100 MB)으로
  올려 두어 긴 메모 셀이 `csv.Error`를 내지 않고, 그래도 나는 `csv.Error`
  (따옴표 불량 등)는 파일명·행 번호를 담은 `InputError`로 바꿉니다.
- `report_row_stats`: 헤더와 셀 수가 다른 행을 `CSV_ROW_LENGTH` 오류로
  (행 번호·셀 수 나열, `LIST_LIMIT = 10`).

### 7.3 `validate.py` -- 1패스 검증 + 쓰기

`process_csv(path, codec, opts, result, out_path=None, row_hook=None) -> details`
가 `validate`와 `clean`의 공용 루프입니다. 행을 한 번만 읽으면서:

1. 열 해석(`resolve_columns`) → `CSV_ID_MISSING` 오류, `CSV_SEQUENCE_MISSING`
   경고.
2. 행마다: 셀 수 불일치면 건너뛰고 통계에 기록(출력에 쓰지 않음); 표준 열의
   결측 집계(`_Agg`), `value`의 숫자 분류, `unit` 집합, `id` 중복, 서열 열의
   알파벳·핵산·소문자 검사(`residues.check_protein_alphabet` 재사용);
   `--drop-missing`이면 표준 열에 결측이 있는 행 제외(미매핑 열의 빈칸은
   무시 -- 메모 열 때문에 행이 지워지지 않게); `row_hook`(F22)이 있으면 호출.
3. `out_path`가 있으면 표준 열(표준 이름, 표준 순서) + 미매핑 열
   (`--drop-unmapped`면 제외)을 UTF-8·LF·`,`·최소 따옴표로 씁니다. 값은 읽은
   그대로.
4. 집계를 finding으로: `CSV_MISSING_VALUES`, `CSV_VALUE_NOT_NUMERIC`,
   `CSV_DECIMAL_COMMA`(소수점 쉼표/천 단위 두 해석 모두 안내), `CSV_UNIT_MIXED`,
   `CSV_ID_DUPLICATE`, `CSV_SEQUENCE_ALPHABET`/`_NUCLEIC`(오류)/`_LOWERCASE`(정보),
   `CSV_ROWS_DROPPED`, `CSV_NO_DATA_ROWS`.

`details`에는 `columns{std: source}`, `unmapped[]`, `missing{col: n}`, `units{}`,
`rows_read`/`rows_written`, 그리고 **`output_written`**(bool)이 들어갑니다. CLI는
출력 파일의 존재 여부가 아니라 이 플래그로만 "출력을 썼다"고 판단합니다
(이전 실행의 낡은 파일이 리포트에 올라가는 사고 방지). 쓰기 중 예외가 나면
반쯤 쓴 출력 파일을 지웁니다.

`inspect_csv`/`format_inspect`는 열별 채움 비율과 예시 값, 자동 매핑
(`Unit -> unit`)을 보여 줍니다.

### 7.4 `crosscheck.py` -- F22: 행이 seed를 가리키는가

`load_reference(path)`는 확장자로 PDB(`.pdb`/`.ent`/`.pdb1` → 체인 ID·서열)와
FASTA(레코드 ID·서열)를 구분해 `Reference(kind, ids, sequences)`를 만듭니다.
참조 FASTA가 UTF-8이 아니면 파일명·위치를 담은 `InputError`("참조를 먼저
`fasta clean`"). `CrossRef(ref)`는 `process_csv`의 `row_hook`으로 들어가 행마다
`id` 일치 **또는** `sequence`(또는 `heavy`+`light` 둘 다) 일치를 검사하고,
`report(result)`가 `CSV_CROSSREF_MATCHED` 정보(by id / by sequence 수)와
`CSV_CROSSREF_UNMATCHED` 경고(미매칭 행 목록, `LIST_LIMIT = 10`), 참조에 비교할
것이 없으면 `CSV_CROSSREF_EMPTY` 경고를 냅니다.

## 8. CLI 계층 `abd_prep/cli.py`

`build_parser()`가 `abd-prep GROUP ACTION ...` 2단계 서브파서를 만듭니다.
`_common_options(parser, with_report, with_encoding)`가 `--json/--strict/
--encoding/--report/--no-report/--report-paths`를 붙이고(텍스트 입력을 읽지
않는 `quickstart`는 `with_encoding=False`), 액션마다
`set_defaults(func=cmd_...)`로 구현 함수를 연결합니다. 모든 구현 함수의
서명은 `cmd_x(args: Namespace, argv: Sequence[str]) -> int`입니다(`argv`는
리포트의 `command` 문자열용). 형식 그룹 밖의 유일한 최상위 명령이
`quickstart`입니다(8.1절).

공통 도우미:

- `_require_file(path)`: 없으면 `InputError`.
- `_default_report_path(output)`: `<출력 stem>.report.json`.
- `_emit(report, args, title)`: `--json`이면 리포트 JSON, 아니면
  `format_table(title)`; `--no-report`가 아니고 `--report`가 있으면 저장.
  `clean`류는 `args.report`가 비어 있을 때 기본 경로를 채운 뒤 호출합니다.
- `_extension_warning`(PDB), `textfile.check_extension`(FASTA/CSV).
- 파일이 여럿인 `validate`는 파일마다 `ValidationResult`를 따로 만들어
  `where`에 파일 이름을 붙인 뒤 합치고, 종료 코드는 최악값을 돌려줍니다.

`main(argv)`는 `argparse`의 `SystemExit`(사용법 오류 → 2), `InputError` → 2,
`UnicodeError` → 2(안전망), `OSError` → 2를 한 곳에서 매핑합니다. 라이브러리
코드에서 종료 코드를 직접 다루지 않는 이유입니다.

옵션 상호 배제·필수 조건은 파서 또는 `cmd_*` 앞부분에서 `InputError`로
처리합니다: `--preset`과 `--target/--heavy/--light/--rename-*` 동시 사용,
`--rename-heavy`만 있고 `--heavy`가 없는 경우, `--against`에 입력 파일 여럿,
`--chain-map`만 있고 `--against`가 없는 경우, 입력과 출력이 같은 파일.

### 8.1 `quickstart.py` -- 간단 모드 (교차 형식 명령)

`abd-prep quickstart`는 ABD 첫 실험에 필요한 파일을 한 폴더에 한 번에
만듭니다. **새 규칙은 하나도 없고** 기존 경로를 순서대로 부르는 조합
계층입니다. `run_quickstart(opts, opener, command) -> QuickstartOutcome`:

```
check_preset(preset)          target 파트 + heavy 파트 필수 → framework: heavy+light=igg, heavy만=nanobody
planned_files / check_overwrite   기존 파일이 있으면 --force 없이는 InputError (쓰기 전에 검사)
1. run_split(SplitOptions.from_preset(..., filenames={"target": "target.pdb", "fv": "antibody_fv.pdb"},
             report_path=""))          fetch → clean → validate; 리포트는 쓰지 않고 findings/inputs/details 만 흡수
2. write_seed()                 inspect_file(fv) → Record(H, L 순서) → write_fasta
                                → validate_records(framework) → crosscheck(records, fv, {})   (F14-F18)
3. write_example_csv()          표준 열 헤더 + seed 한 행(자리 표시 값) → process_csv(row_hook=CrossRef(seed))  (F19-F22)
4. readme_text()                README.txt: 출처(RCSB 개정일·해상도), 파일별 용도, 레시피별 매핑, 주의 사항, 재현 명령
report.details = {preset, framework, split, seed, csv, files};  outputs 5개;  <out>/quickstart.report.json
```

finding은 그룹 접두사를 유지합니다(`target/`, `fv/`는 split이, `seed/`,
`csv/`는 `_prefix()`가 붙임). 그래서 리포트가 네 명령을 손으로 돌린 것과
같은 모양이고, 사용자 가이드의 코드 설명을 그대로 재사용합니다.

설계 결정:

- **고정 파일 이름**(`target.pdb`, `antibody_fv.pdb`, ...). split의
  `<name>_<group>.pdb` 규칙은 한 폴더에 프리셋 여러 개를 돌리는 경우를 위한
  것이고, quickstart는 폴더 하나가 실험 하나이므로 초보자가 그대로 타이핑할
  수 있는 이름을 택했습니다. R5(이름이 정체를 주장하면 안 됨)는 이름에
  항체명을 넣지 않고 출처를 README.txt·리포트에 적는 것으로 지킵니다. 이를
  위해 `SplitOptions.filenames`(그룹 → 파일 이름) 재정의를 추가했고,
  `run_split`은 이름에 경로 구분자가 있으면 거부합니다.
- **CSV는 템플릿.** 유일한 행은 seed 자신이며 `value`는 자리 표시용(`1.0 nM`),
  `notes`에 측정값이 아님을 적습니다. 열은 igg면 `heavy`+`light`, nanobody면
  `sequence`. 숫자를 지어내 결과처럼 보이게 하지 않는다는 §1 원칙과, "CSV
  경로를 끝까지 한 번 돌려 본다"는 목적 사이의 타협입니다. 행이 seed와
  일치하므로 `CSV_CROSSREF_MATCHED 1/1`이 나옵니다.
- **seed FASTA는 쓴 Fv 파일에서 다시 읽어** 만듭니다(원본 아님). 그래서
  FASTA == PDB가 구성상 보장되고, `crosscheck`는 그 사실을 리포트에 남기는
  역할입니다(`FASTA_PDB_IDENTITY` 2건). 헤더는 `pdb inspect --fasta`와 같은
  `roles.fasta_header` 경로를 쓰므로 R7 확정 시 상수 한 곳만 바뀝니다.
- **덮어쓰기 거부가 기본.** 사용자가 고친 seed 파일을 재실행이 지우는 사고를
  막습니다. 검사는 어떤 파일도 쓰기 전에 합니다(`check_overwrite` →
  `os.makedirs`).
- **`--encoding` 없음.** 사용자 텍스트를 읽지 않으므로
  `_common_options(with_encoding=False)`.

테스트는 `tests/test_quickstart.py`: 가짜 `opener`로 1N8Z 미니 픽스처를
"내려받아" 전 과정을 돌리는 단위 테스트, 캐시 디렉터리에 미니 파일을 심어
네트워크 없이 CLI를 끝까지 돌리는 테스트, 그리고 실제 1N8Z가 캐시에 있을
때만 도는 회귀 테스트(seed 서열이 트라스투주맙 VH/VL 상수와 일치).

## 9. report.json 스키마

```json
{
  "tool":    {"name": "abd-prep", "version": "0.1.0", "python": "3.9.25"},
  "command": "abd-prep <group> <action> ...",
  "generated_at": "2026-10-26T19:19:35Z",
  "inputs":  [{"kind": "pdb|fasta|csv", "path": "...", "sha256": "...", "bytes": 478}],
  "outputs": [{"kind": "...", "path": "...", "sha256": "...", "bytes": 471}],
  "findings": [{"level": "info|warning|error", "code": "...", "message": "...", "where": "..."}],
  "details": { }
}
```

`details`는 명령별로 다릅니다.

| 명령 | `details` 내용 |
|---|---|
| `pdb clean` | `atoms_in`, `atoms_out`, `hetatm_removed{resname: n}`, `rename{}`, `chains_out{chain: {residues, atoms, first, last}}`, `ter_records_written`, 변환·제거 목록 |
| `pdb validate` | `files{name: {models, ter_records, has_end, hetatm{}, chains{...sequence, nonstandard}}}`, `preset`, 파일별 `preset_part` |
| `pdb split` | `sources[{pdb_id, path, from_cache, bytes, entry{title, revision_date, initial_release_date, resolution_angstrom, experimental_method, model_count}, entry_error}]`, `parts{target, fv: {selection, out_chain, clean{...}, validation{...}}}` |
| `fasta clean` | `framework`, `records[{id, description, length, role, role_source, lineno}]` |
| `fasta validate` | `files{name: {framework, records[], against{chain_map, chains{role: {offset, identical, overlap, mismatches, pdb_gaps, uncovered, extra_n, extra_c}}}}}` |
| `csv clean/validate` | `columns{std: source}`, `unmapped[]`, `missing{col: n}`, `units{unit: n}`, `rows_read`, `rows_written`, `output_written`, `against{matched, by_id, by_sequence, unmatched}` |

`--against`에 쓴 참조 파일도 `inputs[]`에 들어갑니다(kind `pdb` 또는 `fasta`).
새 필드를 추가할 때는 기존 키를 바꾸지 말고 덧붙이세요 -- 사용자 스크립트가
`findings[].code`와 `inputs[].sha256`에 의존합니다.

## 10. 확장 레시피

### 10.1 프리셋 추가

`abd_prep/pdb/presets.py`의 `PRESETS`에 항목을 하나 더합니다. 잔기 수와
결손은 **실제 파일로 측정**해서 적습니다(`pdb split`을 한 번 돌리고 `inspect`
결과를 보면 됩니다).

```python
"my-antigen-myab": Preset(
    name="my-antigen-myab",
    description="...",                                  # README 프리셋 표에도 한 줄
    parts=(
        Part("target", "7XYZ", "A:30-250"),             # 체인 ID 유지 (R1 미확정)
        Part("light",  "7XYZ", "L:1-108", "L"),
        Part("heavy",  "7XYZ", "H:1-121", "H"),
    ),
    expect_residues={"A": 215, "L": 108, "H": 121},     # 측정값
    expected_gaps=frozenset({("A", 120, 128)}),         # (출력 체인, 앞, 뒤) -- 있을 때만
    notes=("번호 체계·출처·주의점",),
),
```

그다음: `tests/test_split.py`에 `validate_parts`/`output_groups` 수준의
단위 테스트를 추가하고(네트워크 없이 가능), 캐시에 원본이 있으면 실행되는
회귀 테스트를 원하면 `PresetRegressionTest` 패턴을 따릅니다. README 프리셋
표와 사용자 가이드 7.4절에 한 줄을 더합니다. `quickstart --preset`은
`PRESETS`를 그대로 읽으므로 따로 등록할 것은 없지만, target 파트와 heavy
파트가 있어야 받아 줍니다(없으면 `check_preset`이 `InputError`).

### 10.2 변형 잔기 변환 규칙 추가

`abd_prep/common/residues.py`의 변환표에 `ModifiedResidue(parent, renames)`
항목을 더합니다. 삭제할 원자 목록은 적지 않습니다 -- 이름 변경 후 모체의
정규 원자 집합에 없는 원자가 자동으로 제거됩니다. 원소 칸을 바꿔야 하는
이름 변경(MSE의 `SE`→`SD`처럼)은 `element_for_atom`이 처리합니다.
`tests/test_pdb_clean.py`에 픽스처 한 줄(합성 잔기)과 테스트를 추가합니다.

### 10.3 finding code 추가

1. 이름: `GROUP_SUBJECT_STATE` 대문자, 그룹 접두는 `FILE_`/`PDB_`/`FASTA_`/`CSV_`.
   예 `PDB_GAP_CA_DISTANCE`, `CSV_UNIT_MIXED`.
2. 수준 선택: ABD가 거부할 것이 확실하거나 도구가 고칠 수 없는 문제 →
   `error`; 사람이 판단해야 할 점 → `warning`; 도구가 한 일·참고 → `info`.
   "맞을 수도 있는" 것을 오류로 만들지 않습니다(사용자가 `--strict`로 경고를
   실패로 올릴 수 있으므로).
3. 메시지에 **위치와 대처**를 넣습니다(`line 5 '1,2' ... Save the file with
   plain '.' decimals`). 코드는 바꾸지 않고 메시지만 다듬는 것은 자유입니다.
4. README의 Finding codes 표(해당 그룹 절), 사용자 가이드 부록 B(자주 만나는
   것이면)에 추가하고, 그 코드를 만드는 테스트를 하나 둡니다
   (`self.assertIn("NEW_CODE", res.codes())`).

### 10.4 액션 추가

`cli.build_parser()`에서 해당 그룹의 서브파서에 `add_parser`로 액션을 붙이고
`_common_options`를 호출한 뒤 `set_defaults(func=cmd_group_action)`. 구현 함수는
`(args, argv) -> int` 서명으로, 파일 존재 확인 → `Report` 생성 → 라이브러리
호출 → `_emit` → `result.exit_code(args.strict)` 순서를 따릅니다. 출력을 쓰는
액션은 `clean`과 같은 리포트 기본 경로 규칙을 지킵니다. `tests/test_cli*.py`에
`run_cli([...]) -> (code, stdout, stderr)` 도우미로 종단 테스트를 추가합니다.

### 10.5 새 형식 그룹 추가

`abd_prep/<group>/` 패키지를 만들고 `common/`만 import합니다(다른 그룹에
의존하지 않음). 최소 구성은 파서(읽기·정규화·쓰기), 검증(`validate_*`
함수가 `ValidationResult`에 덧붙임), 선검사는 `textfile.precheck` 재사용,
CLI에 `clean/validate/inspect` 세 액션. 패키지 디렉터리 이름이 표준
라이브러리 모듈과 겹치지 않는지 확인하세요(`csvdata`, `ranges`가 그래서 생긴
이름입니다). `pyproject.toml`은 `packages.find include=["abd_prep*"]`라 새
패키지가 자동으로 포함되지만, 빌드 후 wheel 내용을 한 번 확인합니다.

### 10.6 ABD 규약이 확인됐을 때

| 확인된 것 | 고칠 곳 |
|---|---|
| 타겟 체인은 `T`여야 한다 | `presets.py`의 target `Part.out_chain`을 `"T"`로; README 프리셋 표·"not yet confirmed" 절; 사용자 가이드 7.4·12장; `tests/test_split.py` 기대 체인 |
| FASTA 헤더는 `>heavy`/`>light`여야 한다 (또는 설명 불허) | `roles.py`의 `ABD_HEADER_ID`(필요하면 `fasta_header`에서 설명 제거); `tests/test_fasta*.py`의 헤더 기대값; README·가이드 |
| CSV 열 이름이 정해져 있다 | `schema.py`의 `STANDARD_COLUMNS`(순서 포함)와 `REQUIRED_COLUMNS`; `tests/test_csv.py`; README 표 |

세 경우 모두 `--rename-target`, 헤더 수동 편집, `--map`으로 사용자 쪽에서도
우회가 가능하므로 릴리스 사이에도 사용자가 막히지는 않습니다.

## 11. 테스트 전략

- **러너는 `unittest`**(의존성 0 원칙). `python3 -m unittest discover -s tests
  -t .` 또는 `pytest tests/`. 한 파일만: `python3 -m unittest
  tests.test_fasta -v`.
- **네트워크 없음.** `fetch`는 `opener` 주입으로 가짜 응답을 받고, 실제 RCSB
  파일을 쓰는 회귀 테스트는 캐시에 파일이 없으면 `skipTest`.
- **픽스처는 생성기로 만든다.** `tests/make_fixtures.py`(PDB; 캐시 원본에서
  잘라 내고 altloc·MSE·MODEL을 주입, `REMARK`로 출처 표기)와
  `tests/make_seq_fixtures.py`(FASTA/CSV; CRLF·BOM·CP949 바이트를 손으로 만들지
  않기 위해). 커밋된 FASTA/CSV 픽스처는 `FixtureProvenanceTest`와
  `make_seq_fixtures.py --check`로 생성기와 바이트 동일함이 고정됩니다.
  픽스처를 바꾸려면 **생성기를 고치고 다시 생성**합니다. 픽스처 목록과
  출처는 README의 Development 절 표에 있습니다.
- **회귀 테스트 두 종류.** (a) 프로토타입 호환: `clean`이 참조 산출물과
  `ATOM`/`HETATM` 레코드가 바이트 동일한지(`test_regression_prototype.py`,
  `test_split.py::PresetRegressionTest`). (b) 리뷰에서 고친 결함 하나당
  테스트 하나(`test_review_fixes.py`, `test_review_fixes_m25.py`) -- 어떤 입력이
  어떤 잘못된 동작을 냈는지 테스트 이름과 docstring에 적습니다.
- **종단 테스트**는 `cli.main([...])`를 직접 호출해 종료 코드·stdout·생성
  파일·`report.json`을 확인합니다(`tests/helpers.py`의 임시 디렉터리 도우미).
  서브프로세스를 띄우지 않으므로 빠릅니다.
- 테스트 파일 ↔ 모듈 대응: `test_ranges`/`test_pdbfile`/`test_pdb_clean`/
  `test_pdb_validate`/`test_fetch`/`test_split`/`test_cli`(PDB),
  `test_textfile`/`test_fasta`/`test_fasta_roles`/`test_crosscheck`/`test_csv`/
  `test_cli_fasta_csv`(FASTA/CSV).
- 지원 Python 범위(3.8-3.13)에서 최소 두 버전으로 돌려 보세요. 3.8
  호환을 위해 `list[str]` 같은 내장 제네릭 대신 `typing.List`를, `match` 문은
  쓰지 않습니다.

## 12. 릴리스 절차

1. `abd_prep/__init__.py`의 `__version__`을 올립니다(`pyproject.toml`은
   동적으로 읽음). 사용자에게 보이는 동작이 바뀌었으면 README의 해당 절을
   함께 고칩니다.
2. 테스트: 두 Python 버전에서 전체 스위트, `make_seq_fixtures.py --check`.
3. 빌드: `python3 -m build`(또는 `uv build`) → `dist/abd_prep-X.Y.Z.tar.gz`,
   `abd_prep-X.Y.Z-py3-none-any.whl`. wheel에 `abd_prep/*`만 들어 있고
   `__pycache__`가 없는지, sdist에 `tests/`·픽스처·`LICENSE`·`README.md`가
   들어 있는지 확인(`MANIFEST.in`).
4. 새 가상 환경에 wheel을 설치해 `abd-prep --version`, `pdb validate`,
   `fasta validate`, `csv validate`를 픽스처로 한 번씩 실행.
5. **공개 전 스캔**: 소스 트리·sdist·git 전체 이력에서 개인 경로, 계정 ID,
   내부 호스트·시스템 이름, 실측 데이터가 없는지 확인합니다
   (`git rev-list --all | xargs git grep -n -i -e <용어> ...`). 픽스처 서열은
   공개 PDB 항목 유래 또는 명시된 합성 서열이어야 합니다.
6. 커밋·태그(`vX.Y.Z`)·push. PyPI에 올린다면 `twine upload dist/*`.

## 13. 보안 고려 사항

도구는 로컬 파일을 읽고 쓰며, 네트워크는 RCSB 두 호스트에만 HTTPS로
접속합니다. 신뢰 경계는 **사용자 입력 파일**과 **RCSB 응답** 둘입니다.

- **PDB ID → URL·경로**: `normalize_pdb_id`의 정규식을 통과한 4글자만 쓰므로
  경로 이탈·호스트 변조가 불가능합니다. 호스트는 상수, 스킴은 HTTPS만,
  60초 타임아웃, 512 MiB 다운로드 상한, 응답 머리 검사, `.part`→`os.replace`.
- **`--name`·출력 경로**: `derive_name`/`--name`은 경로 구분자를 거부합니다.
  출력 경로 자체는 사용자가 지정하는 것이므로 제한하지 않습니다.
- **입력 파일 파싱**: 고정 칸 슬라이스·`csv.reader`·단순 상태 기계만 쓰고
  `eval`/`pickle`/동적 import는 없습니다. 크기 상한(100 MB 검사, CSV 필드
  상한)과 `MAX_SEQUENCE_LENGTH`로 비정상 입력을 조기에 거부합니다. 디코딩
  오류는 트레이스백 대신 종료 1/2.
- **리포트의 정보 노출**: 경로 기본값 basename. `command`도 같은 규칙으로
  축약됩니다. `--report-paths full`은 사용자가 명시할 때만.
- **쓰기 실패 시**: 부분 출력 파일을 지우고(CSV) `.part`를 남기지 않습니다
  (fetch). 입력과 출력이 같은 경로면 거부합니다.
- 서드파티 코드가 없으므로 공급망 위험은 Python 표준 라이브러리로 한정됩니다.

## 14. 알려진 한계와 후속 과제

| 항목 | 현재 | 후속 방향 |
|---|---|---|
| mmCIF 입력 | 읽지 않음(`fetch`가 mmCIF 전용 항목을 알려만 줌) | 최소 mmCIF 파서(`_atom_site` 루프만) 또는 선택 extra |
| altloc 선택 단위 | 원자별 점유율 | 잔기 전체를 한 컨포머로 맞추는 옵션(`--altloc-by-residue`) |
| Fv 경계 | 프리셋 문헌값, 사용자 지정 | ANARCI(extra)로 PDB·FASTA 공통 자동 산출, 휴리스틱 경고를 실제 절단 제안으로 |
| 프리셋 | 2종(HER2/트라스투주맙) | 다른 공개 항체-항원 복합체 2-3종 |
| `.xlsx` 입력 | 거부 + 안내 | `openpyxl` extra |
| `split` 리포트의 `fv/file` `PDB_HETATM_REMOVED` | 파트 선택 범위가 아닌 **전체 파일** 기준 집계(타겟과 같은 숫자) -- 사실이지만 읽기 어려움 | 파트 선택 범위 안의 HETATM만 세도록 `clean_atoms`에 선택 적용 후 집계 |
| F18 정렬 | 체인당 오프셋 하나(전수 탐색) | 삽입/결실이 있는 다른 항체 비교는 범위 밖으로 유지 |
| 블로그 참조 | README에 링크 없이 언급 | 글이 공개 URL을 가지면 링크 |

## 15. 개발 이력과 주요 결정

이 도구는 두 개의 프로토타입 스크립트 -- `clean_pdb.py`(2-패스 altloc
알고리즘 + 체인 선택 + HETATM 제거 + 체인 재명명, 100줄)와
`inspect_pdb.py`(요약, 45줄) -- 에서 출발했습니다. 블로그의 PyMOL 절차를
스크립트로 옮긴 것이었고 ABD에 실제로 올린 파일 4개(6LBX 타겟, 6B9Z Fv 등)를
만들었습니다. 패키지화하면서 프로토타입의 문제 여섯 가지를 고쳤습니다:
원본 `TER`이 그대로 복사됨(제거한 체인의 `TER`이 남음) → `TER`/`END` 재생성;
삽입 코드 미지원(`int(resSeq)` 비교) → Kabat 번호 대응 범위 문법; HETATM 일괄
삭제가 CSO/MSE 백본 잔기도 지움 → 폴리머형 판정 + 경고 + `--std-residues`;
`MODEL` 미처리 → 첫 모델만; `open()` 누수 → `with`; 빈 선택을 조용히 빈 파일로
씀 → `PDB_SELECTION_EMPTY` 오류. 프로토타입 산출물 4개와의 `ATOM` 레코드
바이트 동일성은 회귀 테스트로 고정했습니다.

단계별 산출물:

1. **PDB 패키지화** -- `common/{residues,result,hashing,report}`,
   `pdb/{pdbfile,ranges,clean,validate}`, `pdb clean|validate|inspect`,
   픽스처 3종. 이때 확인한 사실: 1N8Z 체인 A/B는 삽입 코드 없는 순차 번호이고
   타겟 체인 C의 결손은 580-591 한 곳(CA-CA 6.7 Å); 6B9Z의 CSO는 B175 하나로 Fv
   범위 밖이며 `ACSO`(0.75)/`BCYS`(0.25) 잔기 수준 이질성이라 원자 단위 altloc
   선택만으로는 잔기가 통째로 사라질 수 있어 잔기 이름 선택을 먼저 하도록
   설계를 바꿨습니다.
2. **코드 리뷰 반영** -- UTF-8 오류 위치 계산, 비ASCII 바이트 트레이스백
   (latin-1 왕복으로 해결), 수소 판정 불일치, 재생성 TER 일련번호 충돌 등 8건.
3. **RCSB fetch · split · 프리셋 · `validate --preset`** -- 리포트 파일명을
   `<name>.report.json`으로(같은 디렉터리에서 프리셋 두 개를 돌려도 덮어쓰지
   않게), 프리셋 체인 라벨은 프로토타입 산출물과 동일(Fv `L`/`H`, 타겟 원본
   유지), 항목 해상도 경고 추가.
4. **FASTA·CSV 그룹** -- 서열·값은 바꾸지 않는다는 원칙을 그대로 확장.
   중쇄/경쇄 역할은 추정하지 않음, 끝 `*`만 제거하고 중간 `*`는 오류, F18은
   gap 인식 정렬, CSV `--map`·인코딩·구분자 정규화, `id` 필수·서열 권장.
   나노바디 픽스처는 공개 항목을 받아 오지 않고 1N8Z VH에 낙타류 hallmark
   치환 4개(V37F/G44E/L45R/W47G)를 넣은 **명시된 합성 서열**로, 전장 경고
   픽스처도 합성 리더 + 1N8Z VH + 합성 C-말단으로 구성해 네트워크 없이
   재현되게 했습니다.
5. **두 번째 코드 리뷰 반영** -- 탭 구분 헤더, 다중 끝 `*`, `--replace-ambiguous
   '*=G'` 허용, 큰 CSV 셀 트레이스백, 낡은 출력 파일이 리포트에 오르는 문제,
   `--against`가 아무것도 비교하지 않고 통과하는 문제 등 12건.
6. **문서·라이선스·패키징** -- 영어 README(finding code 전체 표 포함), MIT,
   `pyproject.toml`, sdist/wheel 검증, 공개 전 용어 스캔, 그리고 이 가이드 두
   편.
7. **간단 모드와 튜토리얼** -- `abd-prep quickstart`(8.1절): 첫 ABD 실험에
   필요한 파일을 한 폴더에 한 번에. 기존 경로의 조합이라 finding code는 늘지
   않았고, `SplitOptions.filenames` 재정의 하나가 추가됐습니다. 함께 고친 것:
   `fasta validate --against`의 리포트 `details.against.pdb`가 절대 경로를
   그대로 담아 basename 모드의 경로 비노출 원칙을 깨고 있던 것(basename으로
   수정, 테스트로 고정). 문서는 `docs/easy-tutorial.ko.md`(터미널 초심자용
   20분 경로)와 사용자 가이드 3.0절.

설계 결정 중 다시 검토할 때 맥락이 필요한 것:

- **Biopython을 쓰지 않은 이유.** 초기 설계는 Biopython이었지만 프로토타입이
  이미 표준 라이브러리로 동작했고, 칸 단위 편집이 바이트 동일성을 보장하며,
  설치 없는 실행이 사용자층(연구자)에 더 중요했습니다. `PDBIO`의 altloc
  처리·레코드 손실 문제도 피했습니다.
- **프리셋 기본값이 1N8Z인 이유.** 블로그 절차와 일치하고, 6B9Z 항체는
  1N8Z 트라스투주맙과 프레임워크 여러 자리가 달라(VL `SPILLSASV` vs
  `SPSSLSASV` 등) "트라스투주맙 유래 변이체"로 다뤄야 하기 때문입니다. 두
  번째 프리셋은 실제로 ABD에 올렸던 산출물을 재현하기 위해 남겼습니다.
- **휴리스틱 길이 범위가 구조 검증보다 넓은 이유.** seed 서열은 몇 잔기의
  여분(태그·클로닝 잔기)을 정당하게 가질 수 있지만, 구조 파일의 Fv는 사용자가
  직접 자른 범위라 더 엄격히 볼 수 있습니다.
- **FASTA/CSV `clean`도 오류가 있을 때 출력을 쓰는 이유.** PDB `clean`과의
  일관성, 그리고 사용자가 결과를 보고 판단할 수 있게 하려는 것입니다. 단
  디코딩 자체가 실패하면 출력이 없습니다.

## 16. 코딩 규약

- Python 3.8 호환: `typing`의 `List`/`Dict`/`Optional`, `dataclasses`, f-string.
  `match`, `list[str]`, `str.removeprefix` 등 3.9+ 전용 기능은 쓰지 않습니다.
- 서드파티 import 금지(테스트 포함). 새 의존이 꼭 필요하면 선택 extra로
  분리하고 import 실패 시 안내 메시지를 냅니다.
- docstring·도움말·finding 메시지는 영어(공개 저장소). 모듈 docstring 첫
  부분에 그 모듈이 구현하는 요구사항과 finding code 표를 적어 두는 관행을
  유지합니다.
- 라이브러리 모듈은 `print`하지 않고 `ValidationResult`에 기록합니다. 종료
  코드는 `cli.py`에서만 결정하고, 입력·옵션 오류는 `InputError`로 올립니다.
- 파일을 바꾸는 코드는 "무엇을 바꾸는가"가 한 함수 안에서 보이게 하고,
  바꾸지 않는 것(좌표, 서열 글자, 셀 값)에는 손대지 않습니다. 새 변환은 항상
  옵션 뒤에 둡니다(`--std-residues`, `--replace-ambiguous`처럼 명시적으로).
- 테스트는 결함 하나당 하나, 이름에 증상을 적습니다. 픽스처 바이트는
  생성기로만 바꿉니다.
- 커밋 메시지는 사용자 관점의 변화를 한 줄로(예: `pdb clean: warn when a
  regenerated TER serial collides with an atom serial`).
