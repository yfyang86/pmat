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
| 9 | p-value adjust (BH) | Example: 0.55622667 |


# Install Process

It is the same as `./proc`.
