# Experimental

- BH for p-value adjustment: `padj_bh` in [utils](./src/utils.cpp)
- Add `config.txt`

The format of `config.txt` is:

```
a: 0.60105772
b: 0.
c: -4.0224
d: 0.89
useBartlett: 1
useFoldedNormal: 1
```

While `useFoldedNormal:1`, the output is:

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
| 9 | adjusted p-value via BH(fdr) | Example: 0.55622667 |


# Install Process

It is the same as `./proc`.

# Reference

[BH] Benjamini, Y., and Hochberg, Y. (1995).  Controlling the false discovery rate: a practical and powerful approach to multiple testing.  _Journal of the Royal Statistical Society Series B_, *57*, 289-300.  doi: 10.1111/j.2517-6161.1995.tb02031.x (URL: https://doi.org/10.1111/j.2517-6161.1995.tb02031.x).  <URL:https://www.jstor.org/stable/2346101>.