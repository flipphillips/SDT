# SDT

2.0

![](logo.png)

---

### Signal Detection methods for _Mathemtica._

This package is designed for doing signal detection theory calculations, based on work from:

* N. A. Macmillan & C. D. Creelman (1991). _Detection Theory: A User's Guide,_ Cambridge University Press.
* Smith, J. E. Keith (1982). Simple Algorithms for M-Alternative Forced Choice, _P&P_.


    

### Example use
```
<<SDT`

(* {signal present, response} *)
data = {{True, True}, {False, True} {False, False}, ... };
hf = HitsAndFalseAlarms[data]
    (* {29/45, 2/45} *)

DPrime @@ hf //N
    (* 2.07165 *)

Criterion @@ hf //N
    (* 0.665462 *)
```

There are other methods and options for specifying designs, etc. I'll add more to the `SDT.nb` file soon.

### License
MIT

### More
For more potentially useful stuff, visit [my lab at www.skidmore.edu/~flip](https://www.skidmore.edu/~flip).
