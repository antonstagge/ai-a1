# Grade level E-D

## Question 1

pi =

| P(c1) | P(c2) |
|-------|-------|
| 0.5   | 0.5   |

A =

|      | c1  | c2  |
|------|-----|-----|
| c1   | 0.5 | 0.5 |
| c2   | 0.5 | 0.5 |

B =

|    | h   | t   |
|----|-----|-----|
| c1 | 0.9 | 0.1 |
| c2 | 0.5 | 0.5 |

## Question 2
| 0.5 | 0.5 |
|-----|-----|


## Question 3
| 0.7 |
|-----|
| 0.3 |

## Question 4
This is ok because of the general product rule (chain rule) is used.
The product rule says that:


P(A<sub>1:n</sub>) = P(A<sub>n</sub>|A<sub>1:n-1</sub>) * P(A<sub>1:n-1</sub>)

In this case A<sub>n</sub> = O<sub>t</sub> = o<sub>t</sub>

and A<sub>1:n-1</sub> = X<sub>t</sub>=x<sub>i</sub>,O<sub>1:t-1</sub>

## Question 5

In our implementation we overwrite the delta matrix for each T iteration, therefore the delta matrix will always be a Nx1 matrix.

The delta_index matrix, however, will be a NxT matrix.
