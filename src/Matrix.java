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

    public Matrix getCol(int col) {
        Matrix new_mat = new Matrix(rows, 1);
        for (int i = 0; i < rows; i++) {
            new_mat.set(i,0, this.get(i, col));
        }
        return new_mat;
    }

    public Matrix getRow(int row) {
        Matrix new_mat = new Matrix(1,cols);
        for (int i = 0; i < cols; i++) {
            new_mat.set(0,i, this.get(row, i));
        }
        return new_mat;
    }

    public Matrix transpose() {
        Matrix new_mat = new Matrix(this.cols, this.rows);
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                new_mat.set(j,i, this.get(i, j));
            }
        }
        return new_mat;
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

    public Matrix scalarMult(double scalar) {
        Matrix new_mat = new Matrix(rows, cols);
        for (int i = 0; i < rows*cols; i++) {
            new_mat.mat[i] = mat[i]*scalar;
        }
        return new_mat;
    }

    public Matrix hadamardProd(Matrix other) {
        if (this.rows != other.rows || this.cols != other.cols) {
            throw new IndexOutOfBoundsException("Need to have same shape");
        }

        Matrix new_mat = new Matrix(rows, cols);
        for (int i = 0; i < this.rows; i++) {
            for (int j = 0; j < this.cols; j++) {
                double val = this.get(i,j)*other.get(i,j);
                new_mat.set(i,j, val);
            }
        }

        return new_mat;
    }

    public Matrix appendCol(Matrix colVec) {
        if (colVec.cols != 1 || colVec.rows != this.rows) {
            throw new IndexOutOfBoundsException("appendCol colVec has to be a vector");
        }
        Matrix new_mat = new Matrix(rows, cols+1);
        for (int x = 0; x < rows; x++) {
            for (int y = 0; y < cols; y++) {
                new_mat.set(x,y,this.get(x,y));
            }
        }
        for (int x = 0; x < rows; x++) {
            new_mat.set(x, cols, colVec.get(x, 0));
        }
        return new_mat;
    }

    public Matrix extractMax() {
        Matrix result = new Matrix(rows, 2);
        for (int r = 0; r < rows; r++) {
            int max_idx = 0;
            double max = 0.0;
            for (int c = 0; c < cols; c++) {
                if (this.get(r, c) > max) {
                    max_idx = c;
                    max = this.get(r,c);
                }
            }
            result.set(r, 0, max);
            result.set(r, 1, max_idx);
        }
        return result;
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

    public double diff(Matrix other) {
        if (this.cols != other.cols || this.rows != other.rows) {
            throw new IndexOutOfBoundsException("matrix size no match in diff");
        }
        double sum = 0.0;
        for (int i = 0; i < this.rows*this.cols; i++) {
            sum += Math.sqrt(Math.pow(this.mat[i]-other.mat[i], 2));
        }
        return sum;
    }
}
