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

    public static void main(String[] args) {
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
}
