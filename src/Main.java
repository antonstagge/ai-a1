public class Main {

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

        // alpha_1
        double[] alpha_temp = new double[A.rows];
        for (int i = 0; i < A.rows; i++) {
            alpha_temp[i] = B.get(i, sequence[0])*pi.get(0,i);
        }

        Matrix alpha = new Matrix(alpha_temp.length,1, alpha_temp);

        // alpha_t
        int T = sequence.length;
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
            //System.out.println(delta);
            //System.out.println("New round inner");
            for (int j = 0; j < N; j++) {
                Matrix Bcol = B.getCol(sequence[t]);
                //System.out.println("Bcol");
                //System.out.println(Bcol);
                Matrix Arowtrans = A.getRow(j).transpose();
                //System.out.println("Arowtrans");
                //System.out.println(Arowtrans);
                //System.out.println("deltaJ");
                //System.out.println(delta.get(j,0));
                Matrix tmp = Bcol.hadamardProd(Arowtrans);
                Matrix tmp2 = tmp.scalarMult(delta.get(j,0));
                big = big.appendCol(tmp2);
            }
            //System.out.println(big);
            Matrix max_result = big.extractMax();
            //System.out.println(max_result.getCol(1));
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

    public static void main(String[] args) {
        //HMM0();
        //HMM1();
        HMM2();
    }
}
