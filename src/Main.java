import java.util.Arrays;

public class Main {

    public static double[] normlizer;

    public static Matrix readMatix(Kattio io) {
        int a_rows = io.getInt();
        int a_cols = io.getInt();
        double[] a_m = new double[a_rows*a_cols];
        for (int i = 0; i< a_rows*a_cols; i++) {
            double d = io.getDouble();
            a_m[i] = d;
        }
        Matrix A = new Matrix(a_rows, a_cols, a_m);
        return A;
    }

    public static int[] readRow(Kattio io) {
        int size = io.getInt();
        int[] seq = new int[size];
        for (int i = 0; i < size; i++) {
            seq[i] = io.getInt();
        }
        return seq;
    }

    public static Matrix AlphaPass(Matrix A, Matrix B, Matrix pi, int[] sequence) {

        int T = sequence.length;
        int N = A.rows;
        normlizer = new double[T];
        Matrix alpha = new Matrix(N, 0);

        // alpha_1
        Matrix current_alpha = B.getCol(sequence[0]).hadamardProd(pi.transpose());
        double new_scale = 0;
        for (int i = 0; i < N; i++) {
            new_scale += current_alpha.get(i, 0);
        }
        normlizer[0] = new_scale;
        current_alpha = current_alpha.scalarMult(1.0/new_scale);

        alpha = alpha.appendCol(current_alpha);


        // alpha_t
        for (int t = 1; t < T; t++) {
            Matrix tmp1 = (current_alpha.transpose().mult(A)).transpose();
            Matrix bcol = B.getCol(sequence[t]);
            Matrix tmp2 = tmp1.hadamardProd(bcol);
            current_alpha = tmp2;

            new_scale = 0;
            for (int i = 0; i < N; i++) {
                new_scale += current_alpha.get(i, 0);
            }
            normlizer[t] = new_scale;
            current_alpha = current_alpha.scalarMult(1.0/new_scale);
            alpha = alpha.appendCol(current_alpha);
        }

        return alpha;
    }

    public static Matrix AlphaPass2(Matrix A, Matrix B, Matrix pi, int[] sequence) {
        int T = sequence.length;
        int N = A.rows;
        normlizer = new double[T];
        Matrix alpha = new Matrix(N, T);

        // alpha_1
        for (int i = 0; i <= N-1; i++) {
            alpha.set(i, 0, pi.get(0, i)*B.get(i, sequence[0]));
            normlizer[0] += alpha.get(i, 0);
        }

        for (int i = 0; i < N; i++) {
            alpha.set(i, 0, alpha.get(i, 0)*(1.0/normlizer[0]));
        }

        for (int t = 1; t <= T-1; t++) {
            normlizer[t] = 0;
            for (int i = 0; i <= N-1; i++) {
                alpha.set(i, t, 0);
                for (int j = 0; j <= N-1; j++) {
                    alpha.set(i, t,
                            alpha.get(i, t) + alpha.get(j, t-1)*A.get(j, i)
                    );
                }
                alpha.set(i, t, alpha.get(i, t)*B.get(i, sequence[t]));
                normlizer[t] += alpha.get(i, t);
            }
            for (int i = 0; i <= N-1; i++) {
                alpha.set(i, t, alpha.get(i, t)*(1.0/normlizer[t]));
            }
        }

        return alpha;

    }

    public static Matrix BetaPass(Matrix A, Matrix B, int[] sequence) {

        int T = sequence.length;
        int N = A.rows;
        Matrix beta = new Matrix(N, 0);

        // beta_1
        double[] ones = new double[N];
        Arrays.fill(ones, 1);
        Matrix current_beta = new Matrix(N, 1, ones);
        current_beta = current_beta.scalarMult(1.0/normlizer[0]);

        beta = beta.appendCol(current_beta);

        for (int t = T-2; t > -1; t--) {
            Matrix right_side = current_beta.hadamardProd(B.getCol(sequence[t+1]));
            current_beta = A.mult(right_side);

            //normalize
            current_beta = current_beta.scalarMult(1.0/normlizer[t]);

            beta = beta.appendCol(current_beta);
        }

        return beta;
    }

    public static Matrix BetaPass2(Matrix A, Matrix B, int[] sequence) {
        int T = sequence.length;
        int N = A.rows;

        Matrix beta = new Matrix(N, T);

        // beta_one
        for (int i = 0; i <= N-1; i++) {
            beta.set(i, T-1, 1.0/normlizer[0]);
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
                        beta.get(i, t)*(1.0/normlizer[t])
                );
            }
        }
        return beta;
    }

    public static void HMM0() {
        Kattio io = new Kattio(System.in, System.out);

        Matrix A = readMatix(io);
        Matrix B = readMatix(io);
        Matrix pi = readMatix(io);

        Matrix temp = pi.mult(A);
        Matrix result = temp.mult(B);

        StringBuilder ss = new StringBuilder();
        ss.append(result.rows).append(" ").append(result.cols).append(" ");
        for(int i = 0; i < result.cols*result.rows; i++){
            ss.append(String.format("%.3f",result.mat[i])).append(" ");
        }
        System.out.println(ss.toString().trim());

        io.close();
    }

    public static void HMM1() {
        Kattio io = new Kattio(System.in, System.out);

        Matrix A = readMatix(io);
        Matrix B = readMatix(io);
        Matrix pi = readMatix(io);
        int[] sequence = readRow(io);

        int T = sequence.length;
        int N = A.rows;

        // alpha_1
        Matrix alpha = B.getCol(sequence[0]).hadamardProd(pi.transpose());

        // alpha_t
        for (int t = 1; t < T; t++) {
            Matrix tmp1 = (alpha.transpose().mult(A)).transpose();
            Matrix bcol = B.getCol(sequence[t]);
            Matrix tmp2 = tmp1.hadamardProd(bcol);
            alpha = tmp2;
        }

        // last step
        double sum = 0.0;
        for (int i = 0; i < alpha.rows; i++) {
            sum += alpha.get(i,0);
        }

        System.out.println(sum);

        io.close();
    }

    public static void HMM2() {
        Kattio io = new Kattio(System.in, System.out);

        Matrix A = readMatix(io);
        Matrix B = readMatix(io);
        Matrix pi = readMatix(io);
        int[] sequence = readRow(io);
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

        for (int i = 0; i < T; i++) {
            System.out.print((int)path[i] + " ");
        }
        System.out.println();
        io.close();
    }

    public static void HMM3() {
        Kattio io = new Kattio(System.in, System.out);

        // Initialization
        Matrix A = readMatix(io);
        Matrix B = readMatix(io);
        Matrix pi = readMatix(io);
        int[] sequence = readRow(io);
        int T = sequence.length;
        int N = A.rows;
        int M = B.cols;

        int maxIter = 50;
        int iters = 0;
        double oldLogProb;
        double logProb = -Double.MAX_VALUE;

        // outer convergence loop
        do {
            oldLogProb = logProb;
            //alpha
            Matrix alpha = AlphaPass2(A, B, pi, sequence);
            //beta
            Matrix beta = BetaPass2(A, B, sequence);
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

            // for t = T-1
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
                    A.set(i,j,nominator/denominator);
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
                    B.set(i, j, nominator/denominator);
                }
            }

            // Compute log probability
            logProb = 0;
            for (int t = 0; t <= T-1; t++) {
                logProb += Math.log(1.0/normlizer[t]);
            }
            logProb *= -1;
            iters++;

        } while (iters < maxIter && logProb > oldLogProb); //&& Math.abs(logProb - oldLogProb) > 0.001);

        StringBuilder ss = new StringBuilder();
        ss.append(A.rows).append(" ").append(A.cols).append(" ");
        for(int i = 0; i < A.cols*A.rows; i++){
            ss.append(String.format("%.5f", A.mat[i])).append(" ");
        }
        ss.append("\n");

        ss.append(B.rows).append(" ").append(B.cols).append(" ");
        for(int i = 0; i < B.cols*B.rows; i++){
            ss.append(String.format("%.5f", B.mat[i])).append(" ");
        }
        System.out.println(ss.toString().trim());
    }

    public static void main(String[] args) {
        //HMM0();
        //HMM1();
        //HMM2();
        HMM3();
    }
}
