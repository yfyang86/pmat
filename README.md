# pmat
PMAT detects the differentially methylated regions (DMRs).

| Contents | Logs |
|:---------|:-----|
| Date | 2022-11-22 | 
| Version | v-0.1 |

## C/C++ program

C/C++ version could be found in `./proc/` folder. This projects extends [metilene](https://www.bioinf.uni-leipzig.de/Software/metilene/).

### Pre-install

Required:

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

### INSTALL

Currently, users could use the dev-version, i.e. compile the source to use this tool:

```bash
git clone https://github.com/yfyang86/pmat
cd ./proc/
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

> **Note**: 
For Windows user, one could try to compile and run `pmat` with Cygwin or WSL2.
For MAC OS user, `Home brew` is recommended to install GSL with `brew install gsl` and properly tune the thread numbers. Otherwise, one should compile the GSL from the source and run `pmat` with `-t 1` option. 

# Examples:

** TODO **

## License

C/C++ version (in the proc folder)
- `mm2.h`, `mm2.cpp`： Apache-2.0
- Others： GPL v2.0

## Change logs:

- [x] Initial the projects


## Wishing List

