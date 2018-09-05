import java.util.ArrayList;

public class Matrix {
    public int rows;
    public int cols;
    public ArrayList<Double> mat;

    public Matrix() {
        rows = 0;
        cols = 0;
        mat = new ArrayList<Double>();
    }

    public Matrix(int r, int c, ArrayList<Double> m) {
        rows = r;
        cols = c;
        mat = m;
    }

    public Double get(int x, int y) {
        
    }
}
