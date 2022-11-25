## Linux

Build Env:

> Ubuntu: 18.04
> gcc: 7.5
> gsl: libgsl-dev (2.4+dfsg-6)


Install

```bash
cd ./Linux64-gcc7
tar vxf metilene-0.1-linux64-gcc7.tar.gz

#debian:
sudo apt install libgsl-dev
./metilene --help

#centos/Redhat
sudo yum install libgsl-devel
./metilene --help
```

## OSX-ARM

Build Env:
> System: version 12.5
> gcc/clang: arm64-apple-darwin21.6.0

### pre-install

Install X code or run `xcode-select --install` in any terminal.

### install

```bash
cd ./release/OSX-M1
tar vxf metilene-0.1-osx-arm64-apple-darwin21.6.0.tag.gz 
./metilene --help
```

## OSX-X86

```bash
cd ./release/OSX-Intel
tar vxf metilene-0.1-osx-x86_64-apple-darwin21.6.0.tar.gz
./metilene --help
```

## Windows

Thanks to the R-core team, now we could install 

Preinstall:  [Rtools42](https://cran.r-project.org/bin/windows/Rtools/rtools42/files/rtools42-5355-5357.exe)

Environment: run `C:\rtools42\ucrt64.exe`, and in the bash, run:

```shell
# to where you download pmat
cd pmat
cd ./release/Win64
tar vxf metilene-0.1-win64-rtools42.tar.gz
./metilene --help
```
