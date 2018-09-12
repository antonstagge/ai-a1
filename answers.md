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
If using the absoulte difference between each element in the matrices you get:
diff A : 0.6077394011548561
diff B : 0.5159043660292578
diff pi: 2.0

However, it is not fair to measure the differnce between the matrices one by one, since they all make up the model together. One could say that A is dependent of B and vice versa. Also, the states in the matrices are only represented by indices. By this we mean that row 0 for the original matrix could have been represented by row 1 in our model.

We thought of another way of maesuring the model fitness: to compute the Viterbi algorithm for the two models.
Meaning to compute the most likely sequence of *states* given the sequence of observations for each model. We can them count at how many time-steps these two paths differ to get a measurement of fitness.
However, this will not work for the same reason as above. The numbers the states are given in the algorithm might not match eachother.   

Another way to measure fitness is just to use the likelihood of the model. Meaning the sum of the last alpha-pass vector. However, this is a measurement of how good the model fit the training data. To solve this one could use one set of data to train on and one set of data to test on.

When training our model on the 10000 observation sequence and using max iterations 20000 it converges after 14117 iterations. By testing both this model and the original model on the given sequence that is 1000 oservations, we get the two log probabilities:
New model Log prob : -1347.7165739029786
Original model Log prob : -1342.9301011102039

Which means that our new model is almost as good as the original since higher log probability is better.
ex: e<sup>-2</sup> > e<sup>-10</sup>
## Question 9
With maxiterations set to 20000 and training on the 10000 long sequence and testing on the 1000 long sequence:

Original model Log prob : -1342.930

| # states | Log Prob |
|-----|-----|
| 1 | -1388.751 |
| 2 | -1352.494 |
| 3 | -1347.717 |
| 4 | -1641.464 (did not converge) |
| 5 | -1664.802 (did not converge) |


This experiment shows that having 3 hidden states is the optimal in this case, which is just as
many as there are in the original model. This might not be the case always. If the structure
of the original model is that one column of the A matrix is all 0s (or very close to 0), it means
that it's very unlikely you'll ever transition into that state. In this case a model with only
2 hidden states might more accurately represent the original model.

To determine the optimal number of hidden states one could do use the same approach as we did above:
Train all the models on one set of data and then test them on another and see which one
has the best performance.
If you don't have much data you might want to consider an approach like k-fold cross-validation.
This is where you divide the data into groups, or folds, and train on all folds except one which you use for testing. You then repeat this but change which fold you leave out.

## Question 10


















d
