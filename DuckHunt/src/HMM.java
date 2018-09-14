import java.util.Arrays;

/**
 * Created by anton on 2018-09-14.
 */
public class HMM {
    public Matrix A;
    public Matrix B;
    public Matrix pi;
    private int[] sequence;
    private int N;
    private int M;
    private int T;
    private double[] normalizer;

    int nr_alpha = 0;
    int nr_train = 0;

    HMM() {

    }

    HMM(HMM copy) {
        this.A = new Matrix(copy.A.rows, copy.A.cols, Arrays.copyOf(copy.A.mat, copy.A.mat.length));
        this.B = new Matrix(copy.B.rows, copy.B.cols, Arrays.copyOf(copy.B.mat, copy.B.mat.length));
        this.pi = new Matrix(copy.pi.rows, copy.pi.cols, Arrays.copyOf(copy.pi.mat, copy.pi.mat.length));
        this.N = copy.N;
        this.M = copy.M;
    }

    HMM(int N, int M) {
        A = Matrix.almostUniformMatrix(N, N);
        B = Matrix.almostUniformMatrix(N, M);
        pi = Matrix.almostUniformMatrix(1, N);
        this.N = N;
        this.M = M;
    }

    private Matrix AlphaPass() {
        nr_alpha++;
        normalizer = new double[T];
        Matrix alpha = new Matrix(N, T);

        // alpha_1
        for (int i = 0; i <= N-1; i++) {
            alpha.set(i, 0, pi.get(0, i)*B.get(i, sequence[0]));
            normalizer[0] += alpha.get(i, 0);
        }

        normalizer[0] = 1.0/normalizer[0];

        for (int i = 0; i < N; i++) {
            alpha.set(i, 0, alpha.get(i, 0)*(normalizer[0]));
        }

        for (int t = 1; t <= T-1; t++) {
            normalizer[t] = 0;
            for (int i = 0; i <= N-1; i++) {
                alpha.set(i, t, 0);
                for (int j = 0; j <= N-1; j++) {
                    alpha.set(i, t,
                            alpha.get(i, t) + alpha.get(j, t-1)*A.get(j, i)
                    );
                }
                alpha.set(i, t, alpha.get(i, t)*B.get(i, sequence[t]));
                normalizer[t] += alpha.get(i, t);
            }

            normalizer[t] = 1.0/normalizer[t];

            for (int i = 0; i <= N-1; i++) {
                alpha.set(i, t, alpha.get(i, t)*(normalizer[t]));
            }
        }

        return alpha;
    }

    /*
    public static double[] Viterbi(Matrix A, Matrix B, Matrix pi, int[] sequence) {
        int T = sequence.length;
        int N = A.rows;

        // delta_1
        double[] delta_temp = new double[A.rows];
        for (int i = 0; i < A.rows; i++) {
            delta_temp[i] = B.get(i, sequence[0])*pi.get(0,i);
        }

        Matrix delta = new Matrix(delta_temp.length,1, delta_temp);
        Matrix delta_idx = new Matrix(N, 0);

        // delta_t
        for (int t = 1; t < T; t++) {
            Matrix big = new Matrix(N,0);
            for (int j = 0; j < N; j++) {
                Matrix Bcol = B.getCol(sequence[t]);
                Matrix Arowtrans = A.getRow(j).transpose();
                Matrix tmp = Bcol.hadamardProd(Arowtrans);
                Matrix tmp2 = tmp.scalarMult(delta.get(j,0));
                big = big.appendCol(tmp2);
            }
            Matrix max_result = big.extractMax();
            delta = max_result.getCol(0);
            delta_idx = delta_idx.appendCol(max_result.getCol(1));
        }

        double max = 0;
        int max_idx = 0;
        for (int i = 0; i < N; i++) {
            if (delta.get(i, 0) > max) {
                max = delta.get(i,0);
                max_idx = i;
            }
        }

        double[] path = new double[T];
        path[T-1] = max_idx;
        for (int t = T-2; t >= 0; t--) {
            int prev = (int) path[t+1];
            path[t] = delta_idx.get(prev, t);
        }

        return path;
    }
    */

    private Matrix BetaPass() {

        Matrix beta = new Matrix(N, T);

        // beta_one
        for (int i = 0; i <= N-1; i++) {
            beta.set(i, T-1, normalizer[0]);
        }

        // beta pass
        for (int t = T-2; t >= 0; t--) {
            for (int i = 0; i <= N-1; i++) {
                beta.set(i, t, 0.0);
                for (int j = 0; j <= N-1; j++) {
                    beta.set(i, t,
                            beta.get(i,t) + A.get(i,j)*B.get(j, sequence[t+1])*beta.get(j, t+1)
                    );
                }
                beta.set(i, t,
                        beta.get(i, t)*(normalizer[t])
                );
            }
        }
        return beta;
    }


    public void updateSequence(int[] obsSeq) {
        T = obsSeq.length;
        this.sequence = obsSeq;
    }
    public void trainHMM(int maxIter, int[] obsSequence) {
        nr_train++;
        nr_alpha = 0;
        T = obsSequence.length;
        this.sequence = obsSequence;
        int iters = 0;
        double oldLogProb;
        double logProb = -Double.MAX_VALUE;

        // outer convergence loop
        do {
            oldLogProb = logProb;
            //alpha
            Matrix alpha = this.AlphaPass();
            //beta
            Matrix beta = this.BetaPass();

            // Gamma
            Matrix[] di_gamma = new Matrix[T];
            for (int t = 0; t < T; t++) {
                di_gamma[t] = new Matrix(N, N);
            }

            Matrix gamma = new Matrix(N, T);

            for (int t = 0; t <= T-2; t++) {
                double denominator = 0;
                for (int i = 0; i <= N-1; i++) {
                    for (int j = 0; j <= N-1; j++) {
                        denominator += alpha.get(i,t)*A.get(i, j)*B.get(j,sequence[t+1])*beta.get(j,t+1);
                    }
                }
                for (int i = 0; i <= N-1; i++) {
                    gamma.set(i,t, 0.0);
                    for (int j = 0; j <= N-1; j++) {
                        double nominator = alpha.get(i,t)*A.get(i, j)*B.get(j,sequence[t+1])*beta.get(j,t+1);
                        di_gamma[t].set(i, j, nominator/denominator);
                        double temp = gamma.get(i, t) + di_gamma[t].get(i,j);
                        gamma.set(i, t, temp);
                    }
                }
            }

            // special case for t = T-1
            double denominator = 0;
            for (int i = 0; i <= N-1; i++) {
                denominator += alpha.get(i, T-1);
            }

            for (int i = 0; i <= N-1; i++) {
                gamma.set(i, T-1, alpha.get(i, T-1)/denominator);
            }

            // NEW ESTIMATES:
            // re-estimate pi

            for (int i = 0; i <= N-1; i++) {
                pi.set(0, i, gamma.get(i, 0));
            }

            // re-estimate A
            for (int i = 0; i <= N-1; i++) {
                for (int j = 0; j <= N-1; j++) {
                    double nominator = 0.0;
                    denominator = 0.0;
                    for (int t = 0; t <= T-2; t++) {
                        nominator += di_gamma[t].get(i,j);
                        denominator += gamma.get(i,t);
                    }
                    if (denominator == 0) {
                        if (i == j) {
                            A.set(i,j, 1.0);
                        } else {
                            A.set(i,j, 0.0);
                        }
                    } else {
                        A.set(i,j,nominator/denominator);
                    }
                }
            }

            // re-restimate B
            for (int i = 0; i <= N-1; i++) {
                for (int j = 0; j <= M-1; j++) {
                    denominator = 0.0;
                    double nominator = 0.0;
                    for (int t = 0; t <= T-1; t++) {
                        if (sequence[t] == j) {
                            nominator += gamma.get(i,t);
                        }
                        denominator += gamma.get(i,t);
                    }
                    if (denominator == 0) {
                        B.set(i, j, 1/M);
                    } else {
                        B.set(i, j, nominator/denominator);
                    }
                }
            }

            // Compute log probability
            logProb = 0;
            for (int t = 0; t <= T-1; t++) {
                logProb += Math.log(normalizer[t]);
            }
            logProb *= -1;
            iters++;


        } while (iters < maxIter && logProb > oldLogProb);
    }

    public double calcProbNextObs(int nextObs) {

        Matrix alpha = this.AlphaPass();

        Matrix nextAlphaStep = new Matrix(N, 1);

        for (int i = 0; i <= N-1; i++) {
            nextAlphaStep.set(i, 0, 0);
            for (int j = 0; j <= N-1; j++) {
                nextAlphaStep.set(i, 0,
                        nextAlphaStep.get(i, 0) + alpha.get(j, T-1)*A.get(j, i)
                );
            }
            nextAlphaStep.set(i, 0, nextAlphaStep.get(i, 0) * B.get(i, nextObs));
        }


        double p_new = 0.0;
        for (int i = 0; i < N; ++i) {
            p_new += nextAlphaStep.get(i, 0);
        }

        return p_new;
    }

    public double speciesFitness(int[] sequence) {
        this.sequence = sequence;
        T = sequence.length;
        Matrix alpha = AlphaPass();

        double logProb = 0;
        for (int t = 0; t <= T-1; t++) {
            logProb += Math.log(normalizer[t]);
        }
        logProb *= -1;

        return logProb;
    }
}
