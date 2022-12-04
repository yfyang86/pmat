# PMAT

| Contents | Logs |
|:---------|:-----|
| Date | 2022-12-04 | 
| Version | v-0.2 |

Pairwise methylation assocation test (PMAT) is a computational tool tailored for identifying DMRs between unordered pairs like twins. In this tool, absolute methylation difference between unorderd pairs was maximized in DMRs identification. A folded normal (FN) test based on the framework of likelihood ratio test was proposed recently and implemented to test the methylation difference between unordered pairs in each methylation region.  In order to improve the approximation precision of a FN test when sample size is not large enough, we further established PMAT with Bartlett correction (PMAT-C). In PMAT-C, a folded normal test with Bartlett correction (FN-C) was implemented to test the methylation differences between unordered pairs.

In PMAT with FN test, the asymptotic distribution of a likelihood ratio statistic under the null hypothesis is a mixture of chi-squared distribution. That is

$$LRT \sim 0.5 \chi_0^2 + 0.5 \chi_1^2.$$

In PMAT-C with FN-C test, the mixture proportion of asymptotic distribution with an empirical proportion $\hat{EP}$ is : 
$$LRT \sim (1-\hat{EP}) \chi_0^2 + \hat{EP} \chi_1^2,$$ where $\hat{EP}$ is estimated using $$\hat{EP} = 0.60105772-4.0224 n^{(-0.89)}.$$

Specially, PMAT-C with $\hat{EP} =0.5$ leads to PMAT.

**NOTE:** Currently, all the parameters are hard coded in the `./proc/src/mathematics.cpp` file (the leading 4 entries in the `samParam` array):
```cpp
double foldednomalpvalue(double *a, int m, double *b, int n, bool logit=true, double logittuning=0.00001, bool nonZero = true){
    bool useBartlette = true;
    // config.txt
    double samParam[6]={0.60105772, 0., -4.0224, 0.89, 1., 1.};
	...
}
```
Users could change the parameter and re-compile to get the corresponding pmat-c.

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
apt install libgsl-dev
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

**NOTE**: Binary release for Mac/Windows/Linux(Ubuntu 18.04) users is also [available](https://github.com/yfyang86/pmat/releases). But due to the library/platform dependencies/limitations, some function may not work as designed (eg, parallel computing on Mac OSX). Please contact the maintainer or raise an [issue](https://github.com/yfyang86/pmat/issues)  if you encounter such problems. Currently, only V0.0.1 is available.

### Test

Modify the input file path `../examples/test.dat` in  `test` proc of the [Makefile](./proc/Makefile) to match your settings:

```
test1: 
	./pmat -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05 ../examples/test.dat > test.mr.DMR
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
./pmat -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05  "$INPUT_FILE" > "$OUTPUT_FILE"
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

```bash
./pmat -t 32 -a A -b B -X 8 -Y 8 -m 5 -d 0.05 ../examples/test.dat > test.mr.DMR
```
Or use the binary release. The top 10 lines in the output file are:

```
chr1    661864  661928  0.066676        5       0.0021516       0.75823 0.79141
chr1    662608  662692  0.113871        5       0.33195 0.70271 0.73315
chr1    713375  713453  0.091691        6       0.5     0.67901 0.69787
chr1    842293  842431  0.242725        5       0.5     0.27554 0.32715
chr1    12505   12754   0.119788        8       0.00046499      0.83472 0.81572
chr1    545078  545215  0.371799        11      0.0082285       0.60262 0.71239
chr1    854528  854739  0.159495        5       0.00016186      0.75215 0.67748
chr1    837884  838368  0.106765        14      0.12877 0.84016 0.81546
chr1    848507  848646  0.076560        5       0.5     0.67124 0.65959
chr1    136631  136718  0.140452        6       0.015472        0.85088 0.83067
```
The following table explains the columns in output file.

| Column| Comment | Example |
|:--------|:--------|:--------|
| 1 | chromosome | Example: chr1 |
| 2| start  | Example: 661864 |
| 3 | end | Example: 661928 |
| 4 | avarge absolute 5mC difference between pairs| Example: 0.066676 |
| 5| #CpGs in a region | Example: 5 |
| 6 | p-value via FN-C test or FN test | Example: 0.0021516 |
| 7 | average 5mC in one group| Example: 0.75823 |
| 8 | average 5mC in the other group | Example: 0.79141 |


Run the following command to add column names and FDR value via BH method. (Perl script is available upon request.)

```bash
echo -e "chr\tstart\tstop\tabs.methyl.diff\tCpGs\tpFN\t5mC.A\t5mC.B" > test.mr.DMR.fdr
num_lines.pl test.mr.DMR > o
pvalues_BP_correction.pl o 0 8 | sort -n -k1 | mycut.pl -v -f1 | format_tab.pl >> test.mr.DMR.fdr
```

Alternatively, readers can add column names and FDR values via BH method using the following R script.
```R
data <- read.table("test.mr.DMR", header = F)
colnames(data) <- c("chr", "start", "stop", "abs.methyl.diff", "CpGs", "pFN", "5mC.A", "5mC.B")
data$FDR <- p.adjust(data$pFN, method = "BH")
write.table(data, "test.mr.DMR.fdr", row.names= F, sep = "\t", quote = F)
```



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
- [ ] Config file support for the Bartlett correction : In the comming Version
  - [x] Documentation and Examples.
  - [ ] `Config.txt` support without re-compile. 
- [x] Adjusted p-value support : In the comming Version

