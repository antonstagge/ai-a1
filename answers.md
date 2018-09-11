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

## Question 6

This is because we use the product rule but backwards:

P(A|B) = P(A,B) / P(B)

In our case:

A = X<sub>t</sub>=x<sub>i</sub>, X<sub>t+1</sub>=x<sub>j</sub>

and

B = O<sub>1:T</sub> = o<sub>1:T</sub>

which is exactly what we get when summing the last alpha vector.
The probability of observing the sequence given the model.

# HMM C
## Question 7
The algorithm converges on 1000 observations, but not on 10000 as it doesn't stop within reasonable time (max iteration 10000).

Convergence can be defined as: The algorithm stops in a maximum. This, however, does not *have* to be a global maximum, it can be a local. You have to decide what level of convergence you want, and you can do this by changing the number of max iterations for the baum-welch algorithm.

## Question 8

When using 7000 observations we got:

A =

| 0.713 | 0.025 | 0.263 |
| 0.131 | 0.684 | 0.185 |
| 0.116 | 0.222 | 0.662 |

B =

| 0.701 | 0.187 | 0.096 | 0.017 |
| 0.102 | 0.465 | 0.310 | 0.123 |
| 0.057 | 0.197 | 0.209 | 0.537 |

Pi =

| 0.0 | 1.0 | 0.0 |

which pretty close for the A and B matrices but pi is totally off.


## Question 9




















d
