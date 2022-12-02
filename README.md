# PMAT
PMAT detects the differentially methylated regions (DMRs) for unordered pairs.

| Contents | Logs |
|:---------|:-----|
| Date | 2022-11-22 | 
| Version | v-0.1 |

## C/C++ program

C/C++ version could be found in `./proc/` folder. 

1. **Abilities**
 - `pmat` utilizes a totally different loss-function to tackle the non-order problem in `twin` data, which fails to distinguish `Effect 1 - Effect 2` versus `Effect 2 - Effect 1`;
 - `pmat` provides several optimization routines with different initial value schemes to solve the corresponding MLE to balance the computation accuracy, speed and scalbility (EM / BFGS);
 - `pmat` detects differentially methylated regions between unordered pairs.

2. **Transferability and Extendibility**
 - `pmat` shares the same options as [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/) .
 - In the future, `pmat` plans to provide options to call the original [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/).

3. **Soft Engineering on** [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/)   
 - Bug fix
 - Port `metilene` from ANSI-C to Modern C++11 (partly)

### Pre-install

#### Required
A modern C++ compiler with a GNU library is needed:

```
C++: C++11
Libraries: GSL-2.x
```

- To enable the parallel computing ability (`-t` option), `pthread` support is **required**.

- The `C++11` featured compiler could be found [here](https://en.cppreference.com/w/cpp/compiler_support).
  
- Install `dev`/`devel` version of GSL. For example:

```bash
# debian/ubuntu 
apt-install libgsl-dev
# Centos/RH
yum install libgsl-devel
```

#### Optional

> **Note**: for Windows user, one should intall [Rtools42](https://cran.r-project.org/bin/windows/Rtools/rtools42/files/rtools42-5355-5357.exe) and run `pmat` in the  `C:\rtools42\ucrt64.exe` environment. In the `ucrt64` bash, run
```bash
# to where you download pmat
cd pmat
cd ./release/Win64
tar vxf PMAT-0.1-win64-rtools42.tar.gz
./PMAT --help
```

> **Note**: For MAC OS user, `Home brew` is recommended to install GSL with `brew install gsl` and properly tune the thread numbers. Otherwise, one should compile the GSL from the source and run `pmat` with `-t 1` option. 

### INSTALL

Currently, users could use the development version, i.e. compile `pmat` from source after the `pre-install` procedure:

```bash
git clone https://github.com/yfyang86/pmat
cd ./pmat/proc/
make
```


**NOTE**: Tested Platform: 
- [x] X86/64(Linux/MAC OSX/Windows);
- [x] ARM(MAC M1-ARM64/Linux-ARMv7/ARMv8);
- [x] MIPS/LoongArch(龙芯3A5000, UOS).

**NOTE**: Binary release for Mac/Windows/Linux(Ubuntu 18.04) users is also [available](https://github.com/yfyang86/pmat/releases). But due to the library/platform dependencies/limitations, some function may not work as designed (eg, parallel computing on Mac OSX). Please contact the maintainer or raise an [issue](https://github.com/yfyang86/pmat/issues)  if you encounter such problems.

### Test

Modify the input file path `../examples/test.dat` in  `test` proc of the [Makefile](./proc/Makefile) to match your settings:

```
test1: 
	./PMAT -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05 ../examples/test.dat > test.mr.DMR
```

And run:

```bash
make test
```



# Usage and Examples

The parameters are the same as [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/). For example, `-t 32` means 32 threads are used; `-a A` represents to one group named `A`; `-b B` represents to the other group named `B`; `-X 8` means more than 8 non-missing values are required to estimate missing values in group A; `-Y 8` means more than 8 non-missing values are required to estimate missing values in group B; `-m 5` means more than 5 CpGs are required in a DMR.
A typical example is:

```bash
INPUT_FILE="your_input_file_path"
OUTPUT_FILE="output_file_path" 
./PMAT -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05  "$INPUT_FILE" > "$OUTPUT_FILE"
```

Here `INPUT_FILE` are the input file. The format should be a TAB separated file:

| Columns | Type | Comment |
|:--------|:--------|:--------|
| chr | character | Example: chr1 |
| pos | integer | Example: 10497 |
| A | double | Example: 0.346938775510204 |
| B | double | Example: 0.369565217391304 |

Here `A` and `B` appear in pair. If there are 10 groups, then there should be 10 (`A` and `B`) pairs. The file should be (first three lines):

```
chr	pos	A	B	A	B	A	B	A	B	A	B	A	B	A	B	A	B	A	B	A	B
chr1	10497	0.346938775510204	0.369565217391304	0.362162162162162	0.322727272727273	0.352941176470588	0.486486486486487	0.430769230769231	0.440860215053763	0.439024390243902	0.410526315789474	0.343137254901961	0.40625	0.25974025974026	0.256756756756757	0.296	0.4375	0.377049180327869	0.421487603305785	0.450549450549451	0.385542168674699
chr1	10525	0.92	0.9375	0.940540540540541	0.969298245614035	1	0.972972972972973	0.984615384615385	0.956989247311828	0.963414634146341	0.958333333333333	1	0.958333333333333	0.987012987012987	0.945945945945946	0.944	0.9375	0.975409836065574	0.933884297520661	0.978021978021978	0.975903614457831
chr1	10542	0.96	0.90625	0.978378378378378	0.969298245614035	0.882352941176471	0.945945945945946	0.96875	0.978494623655914	0.963414634146341	0.96875	0.941176470588235	0.957894736842105	0.922077922077922	0.959459459459459	0.96	0.958762886597938	0.959016393442623	0.983471074380165	0.934065934065934	0.951807228915663
```


Run the test code in linux

```
	./PMAT -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05 ../examples/test.dat > test.mr.DMR
```
Or use the binary release. The top 10 lines in the output file are:

```
chr1	12505	12679	1	0.108321	7	0.0034718	0.90758	0.83462	0.80487
chr1	661864	661928	1	0.066676	5	0.0016	0.65065	0.75823	0.79141
chr1	662608	662692	1	0.113871	5	0.28328	0.79049	0.70271	0.73315
chr1	715040	715121	1	0.214549	5	0.46519	0.17554	0.56345	0.69286
chr1	545153	545215	1	0.362152	5	0.032756	0.048372	0.59842	0.76797
chr1	713375	713449	1	0.083285	5	0.46782	0.79534	0.68649	0.72001
chr1	842293	842431	1	0.242725	5	0.41702	0.72288	0.27554	0.32715
chr1	136631	136718	1	0.140452	6	0.014309	0.28271	0.85088	0.83067
chr1	136814	136876	1	0.121317	5	0.45345	0.12209	0.87179	0.82835
chr1	136912	137169	1	0.183106	11	0.33916	0.88037	0.62192	0.65789
```

| Column| Comment | Example |
|:--------|:--------|:--------|
| 1 | chromosome | Example: chr1 |
| 2| start  | Example: 10497 |
| 3 | end | Example: 12679 |
| 4 | q-value | Example: 1 |
| 5 | avarge absolute 5mC difference between pairs| Example: 0.108321 |
| 6| #CpGs in a region | Example: 7 |
| 7 | p-value via FN-C test or FN test | Example: 0.0034718 |
| 8 | p-value via 2d KS test | Example: 0.90758 |
| 9 | average 5mC in one group| Example: 0.83462 |
| 10 | average 5mC in the other group | Example: 0.80487 |


Run the following command to add column names and FDR value via BH method.

```
echo -e "chr\tstart\tstop\tq-value\tabs.methyl.diff\tCpGs\tpFN\tp2DKS\tpre\tpost\tpFN.fdr" > test.mr.DMR.fdr
num_lines.pl test.mr.DMR > o
pvalues_BP_correction.pl o 0 8 | sort -n -k1 | mycut.pl -v -f1 | format_tab.pl >> test.mr.DMR.fdr
```

#

**TODO**
**TODO**
**TODO**

## License

C/C++ version (in the `./proc` folder):
- `mm2.h`, `mm2.cpp`： Apache-2.0
- Others： GPL v2.0

## Change logs:

- [x] 2022-11-22: Initial the pmat projects
- [ ] 2022-11: v-0.1 binary-preview-release 
  - [x] Linux (X86-64)
  - [x] Mac OSX (ARM/Intel)
  - [x] Windows
- [ ] 2022-11: Example/Simulation
	- [x] Simple demo
	- [ ] Full example (in plan)


## Wishing List

- [ ] [Simulation Demos](./TODO.md): 
