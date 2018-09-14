import java.util.Random;

public class Matrix {
    public int rows;
    public int cols;
    public double[] mat;

    public Matrix(int r, int c) {
        rows = r;
        cols = c;
        mat = new double[r*c];
    }

    public Matrix(int r, int c, double[] m) {
        rows = r;
        cols = c;
        mat = m;
    }

    public Double get(int x, int y) {
        if (x > rows-1 || x < 0 || y > cols-1 || y < 0) {
            throw new IndexOutOfBoundsException("Not a valid index in matrix");
        }
        int index = x*cols + y;
        return mat[index];
    }


    public void set(int x, int y, double val) {
        if (x > rows-1 || x < 0 || y > cols-1 || y < 0) {
            throw new IndexOutOfBoundsException("Not a valid index in matrix");
        }
        int index = x*cols + y;
        mat[index]= val;
    }

    public String toString() {
        StringBuilder ss = new StringBuilder();
        for (int x = 0; x < rows; x++) {
            for (int y = 0; y < cols; y++) {
                ss.append(this.get(x, y).toString()).append(", ");
            }
            ss.append("\n");
        }
        return  ss.toString();
    }

    public static Matrix almostUniformMatrix(int N, int M) {
        double[] matrix = new double[N*M];
        Random random = new Random();
        for(int row = 0; row < N; row++){
            double rowSum = 0.0;
            for(int column = 0; column < M; column++){
                if(column == M-1){
                    //Last in row
                    double newValue = 1 - rowSum;
                    matrix[row*M + column] = newValue;
                } else {
                    double temp = (random.nextDouble()-0.5)/100;
                    double newValue = 1.0/M + temp;
                    rowSum += newValue;
                    matrix[row*M + column] = newValue;
                }
            }
        }
        return new Matrix(N, M, matrix);
    }
}
