import java.util.ArrayList;

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

    public Matrix mult(Matrix other) {
        if (cols != other.rows) {
            throw new IndexOutOfBoundsException("Matrices A amount of columns have to match matrix Bs amount of rows to be multiplied");
        }
        int new_cols = other.cols;
        int new_rows = rows;
        Matrix new_mat = new Matrix(new_rows, new_cols);

        for (int x = 0; x < rows; x++) {
            for (int y = 0; y < new_cols; y++) {
                double value = 0.0;
                for (int m = 0; m < cols; m++) {
                    value += this.get(x, m)*other.get(m, y);
                }
                new_mat.set(x, y, value);
            }
        }
        return new_mat;
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
}
