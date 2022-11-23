# pmat
PMAT detects the differentially methylated regions (DMRs).

| Contents | Logs |
|:---------|:-----|
| Date | 2022-11-22 | 
| Version | v-0.1 |

## C/C++ program

C/C++ version could be found in `./proc/` folder. 

1. **Abilities**
 - `pmat` utilizes a totally different loss-function to tackle the non-order problem in `twin` data, which fails to distinguish `Effect 1 - Effect 2` versus `Effect 2 - Effect 1`;
 - `pmat` provides several optimization routines with different initial value schemes to solve the corresponding MLE to balance the computation accuracy, speed and scalbility (EM / BFGS).

2. **Transferability and Extendibility**
 - `pmat` shares the same options as [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/) .
 - In the future, `pmat` plans to provide options to call the original [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/).

3. **Soft Engineering on** [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/)   
 - Bug fix
 - Port `metilene` from ANSI-C to Modern C++11 (partly)

### Pre-install

#### Required
A modern C++ compiler with a GNU library is need:

```
C++: C++11
Libraries: GSL-2.x
```

- To enable the parallel computing ability (`-t` option), `pthread` support is required. 
  
- Install `dev`/`devel` version of GSL. For example:

```bash
# debian/ubuntu 
apt-install libgsl-dev
# Centos/RH
yum install libgsl-devel
```

#### Optional

> **Note**: 
For Windows user, one could try to compile and run `pmat` with Cygwin or WSL2.
For MAC OS user, `Home brew` is recommended to install GSL with `brew install gsl` and properly tune the thread numbers. Otherwise, one should compile the GSL from the source and run `pmat` with `-t 1` option. 

### INSTALL

Currently, users could use the dev-version, i.e. compile the source to use this tool:

```bash
git clone https://github.com/yfyang86/pmat
cd ./pmat/proc/
make
```

Tested Platform: 
- X86/64(Linux/MAC OSX);
- ARM(MAC M1/ARM 64);
- LoongArch(龙芯3A5000, UOS).

### Test

Modify the test proc in the Makefile to match your settings:

```
test1: 
	./metilene -t 32 -a A -b B -X 104 -Y 104 -m 5 -d 0.01 ../../../data/twin.test.data.txt > test.mr.DMR
```

And run:

```bash
make test
```

The parameters are the same as [metline](https://www.bioinf.uni-leipzig.de/Software/metilene/). For example, `-t 32` means 32 threads are used.



# Usage and Examples

**TODO**
**TODO**
**TODO**

## License

C/C++ version (in the `./proc` folder):
- `mm2.h`, `mm2.cpp`： Apache-2.0
- Others： GPL v2.0

## Change logs:

- [x] 2022-11-22: Initial the pmat projects
- [ ] 2022-11: v-0.1 binary-preview-release (in plan)
- [ ] 2022-11: Example/Simulation (in plan)


## Wishing List

- [ ] [Simulation Demos](./TODO.md): 
