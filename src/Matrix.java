
import java.util.*;

public class Matrix {

    private static final double EPS = 1e-9;

    /* =====================================================
       =============== INTERNAL UTILITIES ==================
       ===================================================== */

    /** Deep-copy a 2-D list so the original is never mutated. */
    public static ArrayList<ArrayList<Double>> deepCopy(ArrayList<ArrayList<Double>> m) {
        ArrayList<ArrayList<Double>> copy = new ArrayList<>();
        for (ArrayList<Double> row : m)
            copy.add(new ArrayList<>(row));
        return copy;
    }

    /** Return the n×n identity matrix. */
    public static ArrayList<ArrayList<Double>> identity(int n) {
        ArrayList<ArrayList<Double>> I = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < n; j++)
                row.add(i == j ? 1.0 : 0.0);
            I.add(row);
        }
        return I;
    }

    /** Swap two rows in-place. */
    private static void swapRows(ArrayList<ArrayList<Double>> m, int r1, int r2) {
        ArrayList<Double> tmp = m.get(r1);
        m.set(r1, m.get(r2));
        m.set(r2, tmp);
    }

    /**
     * Subtract (factor * base) from target, element-wise.
     * Used by Gaussian elimination steps.
     */
    private static void eliminate(ArrayList<Double> base,
                                  ArrayList<Double> target,
                                  double factor) {
        for (int i = 0; i < target.size(); i++)
            target.set(i, target.get(i) - factor * base.get(i));
    }

    /* =====================================================
       ================= BASIC OPERATIONS ==================
       ===================================================== */

    /** Element-wise addition of two matrices of the same size. */
    public static ArrayList<ArrayList<Double>> add(
            ArrayList<ArrayList<Double>> A,
            ArrayList<ArrayList<Double>> B) {

        var result = deepCopy(A);
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < A.get(0).size(); j++)
                result.get(i).set(j, A.get(i).get(j) + B.get(i).get(j));
        return result;
    }

    /** Element-wise subtraction (A - B). */
    public static ArrayList<ArrayList<Double>> subtract(
            ArrayList<ArrayList<Double>> A,
            ArrayList<ArrayList<Double>> B) {

        var result = deepCopy(A);
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < A.get(0).size(); j++)
                result.get(i).set(j, A.get(i).get(j) - B.get(i).get(j));
        return result;
    }

    /**Returns true when A and B can be multiplied (A.cols == B.rows)*/
    public static boolean canMultiply(ArrayList<ArrayList<Double>> A,
                                      ArrayList<ArrayList<Double>> B) {
        return A.get(0).size() == B.size();
    }

    /** Matrix multiplication A × B. */
    public static ArrayList<ArrayList<Double>> multiply(
            ArrayList<ArrayList<Double>> A,
            ArrayList<ArrayList<Double>> B) {

        if (!canMultiply(A, B))
            throw new IllegalArgumentException(
                    "Cannot multiply: A has " + A.get(0).size() +
                            " columns but B has " + B.size() + " rows.");

        int rows  = A.size();
        int cols  = B.get(0).size();
        int inner = A.get(0).size();

        ArrayList<ArrayList<Double>> result = new ArrayList<>();
        for (int i = 0; i < rows; i++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int j = 0; j < cols; j++) {
                double sum = 0;
                for (int k = 0; k < inner; k++)
                    sum += A.get(i).get(k) * B.get(k).get(j);
                row.add(sum);
            }
            result.add(row);
        }
        return result;
    }

    /** Transpose of a matrix (rows become columns). */
    public static ArrayList<ArrayList<Double>> transpose(
            ArrayList<ArrayList<Double>> A) {

        int rows = A.size();
        int cols = A.get(0).size();
        ArrayList<ArrayList<Double>> T = new ArrayList<>();
        for (int j = 0; j < cols; j++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int i = 0; i < rows; i++)
                row.add(A.get(i).get(j));
            T.add(row);
        }
        return T;
    }

    /** Scale every element of a matrix by a scalar. */
    public static ArrayList<ArrayList<Double>> scale(
            ArrayList<ArrayList<Double>> A, double scalar) {

        var result = deepCopy(A);
        for (int i = 0; i < A.size(); i++)
            for (int j = 0; j < A.get(0).size(); j++)
                result.get(i).set(j, A.get(i).get(j) * scalar);
        return result;
    }

    /**RREF HELPERS

     Divide every element of a row by the pivot value, making the pivot = 1. */
    private static void normalizeRow(ArrayList<ArrayList<Double>> A,
                                     int pivotRow,
                                     int col,
                                     int cols) {
        double pivot = A.get(pivotRow).get(col);
        for (int j = 0; j < cols; j++)
            A.get(pivotRow).set(j, A.get(pivotRow).get(j) / pivot);
    }

    /**
     * Find the first row at or below startRow where column col is non-zero.
     * Returns -1 if no such row exists (whole column is zero below startRow).
     */
    private static int findPivotRow(ArrayList<ArrayList<Double>> A,
                                    int startRow,
                                    int col,
                                    int rows) {
        for (int i = startRow; i < rows; i++)
            if (Math.abs(A.get(i).get(col)) > EPS)
                return i;
        return -1;
    }

    /**
     * Zero-out every entry in column col EXCEPT the pivot row,
     * using row-addition elimination.
     */
    private static void eliminateColumn(ArrayList<ArrayList<Double>> A,
                                        int pivotRow,
                                        int col,
                                        int rows) {
        for (int i = 0; i < rows; i++) {
            if (i == pivotRow) continue;
            double factor = A.get(i).get(col);
            if (Math.abs(factor) > EPS)
                eliminate(A.get(pivotRow), A.get(i), factor);
        }
    }

    /**RREF**/

    /**
     * Reduced Row Echelon Form (RREF) via Gauss–Jordan elimination.
     * The input matrix is not modified.
     */
    public static ArrayList<ArrayList<Double>> rref(
            ArrayList<ArrayList<Double>> matrix) {

        var A = deepCopy(matrix);
        int rows = A.size();
        int cols = A.get(0).size();
        int pivotRow = 0;

        for (int col = 0; col < cols && pivotRow < rows; col++) {
            int pivotIndex = findPivotRow(A, pivotRow, col, rows);
            if (pivotIndex == -1) continue;
            swapRows(A, pivotRow, pivotIndex);
            normalizeRow(A, pivotRow, col, cols);
            eliminateColumn(A, pivotRow, col, rows);
            pivotRow++;
        }
        return A;
    }

    /**ROW ECHELON (REF)**/

    /**
     * Row Echelon Form (upper-triangular, forward elimination only).
     * zeroed out rows ABOVE the pivot too — that produces RREF, not REF.
     * Here we only eliminate rows BELOW the pivot.
     */
    public static ArrayList<ArrayList<Double>> ref(
            ArrayList<ArrayList<Double>> matrix) {

        var A = deepCopy(matrix);
        int rows = A.size();
        int cols = A.get(0).size();
        int pivotRow = 0;

        for (int col = 0; col < cols && pivotRow < rows; col++) {
            int pivotIndex = findPivotRow(A, pivotRow, col, rows);
            if (pivotIndex == -1) continue;
            swapRows(A, pivotRow, pivotIndex);

            for (int i = pivotRow + 1; i < rows; i++) {
                double factor = A.get(i).get(col) / A.get(pivotRow).get(col);
                if (Math.abs(factor) > EPS)
                    eliminate(A.get(pivotRow), A.get(i), factor);
            }
            pivotRow++;
        }
        return A;
    }

    /**SYSTEM SOLVER **/

    /** Append column vector b to matrix A, creating the augmented matrix [A|b]. */
    public static ArrayList<ArrayList<Double>> augmented(
            ArrayList<ArrayList<Double>> A,
            ArrayList<Double> b) {

        var aug = deepCopy(A);
        for (int i = 0; i < b.size(); i++)
            aug.get(i).add(b.get(i));
        return aug;
    }

    /**
     * Solve the linear system A·x = b using RREF.
     *
     * Returns a map of variable names → values (e.g. "X1", "X2", …).
     * Throws IllegalStateException for no solution.
     * Throws IllegalStateException for infinitely many solutions (free variables).
     */
    public static Map<String, Double> solve(
            ArrayList<ArrayList<Double>> A,
            ArrayList<Double> b) {

        int n = A.get(0).size(); // number of unknowns
        var aug  = augmented(A, b);
        var rref = rref(aug);
        int rows = rref.size();
        int cols = rref.get(0).size(); // n + 1

        // Check for inconsistency or free variables
        for (int i = 0; i < rows; i++) {
            boolean allZeroLHS = true;
            for (int j = 0; j < n; j++) {
                if (Math.abs(rref.get(i).get(j)) > EPS) {
                    allZeroLHS = false;
                    break;
                }
            }
            double rhs = rref.get(i).get(cols - 1);

            if (allZeroLHS && Math.abs(rhs) > EPS)
                throw new IllegalStateException("No solution (inconsistent system).");

            if (allZeroLHS) continue; // zero row — skip

            // Count leading variables in this row
            int pivotCol = -1;
            for (int j = 0; j < n; j++) {
                if (Math.abs(rref.get(i).get(j)) > EPS) { pivotCol = j; break; }
            }
            // If a row has a pivot but non-zero entries in other columns → free variables
            for (int j = pivotCol + 1; j < n; j++) {
                if (Math.abs(rref.get(i).get(j)) > EPS)
                    throw new IllegalStateException(
                            "Infinitely many solutions (free variable: X" + (j + 1) + ").");
            }
        }

        Map<String, Double> result = new LinkedHashMap<>();
        for (int i = 0; i < n; i++)
            result.put("X" + (i + 1), rref.get(i).get(cols - 1));
        return result;
    }

    /* =====================================================
       ================== LU DECOMPOSITION =================
       ===================================================== */

    public static class LUResult {
        public final ArrayList<ArrayList<Double>> L;
        public final ArrayList<ArrayList<Double>> U;

        public LUResult(ArrayList<ArrayList<Double>> L,
                        ArrayList<ArrayList<Double>> U) {
            this.L = L;
            this.U = U;
        }
    }

    /**
     * LU decomposition without pivoting.
     * A = L · U  where L is lower-unit-triangular and U is upper-triangular.
     * Throws if a zero pivot is encountered (use LDU for a pivot check).
     */
    public static LUResult decomposeLU(ArrayList<ArrayList<Double>> matrix) {
        int n = matrix.size();
        var L = identity(n);
        var U = deepCopy(matrix);

        for (int i = 0; i < n; i++) {
            if (Math.abs(U.get(i).get(i)) < EPS)
                throw new IllegalStateException(
                        "Zero pivot at row " + i + " — LU decomposition requires non-zero pivots.");

            for (int k = i + 1; k < n; k++) {
                double factor = U.get(k).get(i) / U.get(i).get(i);
                L.get(k).set(i, factor);
                eliminate(U.get(i), U.get(k), factor);
            }
        }
        return new LUResult(L, U);
    }

    /**LDU DECOMPOSITION **/

    public static class LDUResult {
        public final ArrayList<ArrayList<Double>> L;
        public final ArrayList<ArrayList<Double>> D;
        public final ArrayList<ArrayList<Double>> U;

        public LDUResult(ArrayList<ArrayList<Double>> L,
                         ArrayList<ArrayList<Double>> D,
                         ArrayList<ArrayList<Double>> U) {
            this.L = L;
            this.D = D;
            this.U = U;
        }
    }

    /**
     * LDU decomposition: A = L · D · U
     *   L — lower unit-triangular
     *   D — diagonal matrix of pivots
     *   U — upper unit-triangular
     */
    public static LDUResult decomposeLDU(ArrayList<ArrayList<Double>> matrix) {
        int n = matrix.size();
        var lu = decomposeLU(matrix); // reuse LU logic
        var L  = lu.L;
        var rawU = lu.U;
        var D  = identity(n);
        var U  = deepCopy(rawU);

        for (int i = 0; i < n; i++) {
            double diag = rawU.get(i).get(i);
            if (Math.abs(diag) < EPS)
                throw new IllegalStateException("Zero pivot at row " + i);
            D.get(i).set(i, diag);
            for (int j = i; j < n; j++)
                U.get(i).set(j, rawU.get(i).get(j) / diag);
        }
        return new LDUResult(L, D, U);
    }

    /**DETERMINANT**/

    /**
     * Determinant via LDU decomposition.
     * det(A) = product of all diagonal entries in D.
     */
    public static double determinant(ArrayList<ArrayList<Double>> matrix) {
        if (matrix.size() != matrix.get(0).size())
            throw new IllegalArgumentException("Determinant requires a square matrix.");

        var ldu = decomposeLDU(matrix);
        double det = 1.0;
        for (int i = 0; i < ldu.D.size(); i++)
            det *= ldu.D.get(i).get(i);
        return det;
    }

    /**INVERSE**/

    /**
     * Matrix inverse via Gauss–Jordan on [A | I].
     * Returns A-1, or throws if A is singular.
     */
    public static ArrayList<ArrayList<Double>> inverse(
            ArrayList<ArrayList<Double>> matrix) {

        int n = matrix.size();
        if (n != matrix.get(0).size())
            throw new IllegalArgumentException("Inverse requires a square matrix.");

        // Build [A | I]
        var aug = deepCopy(matrix);
        var eye = identity(n);
        for (int i = 0; i < n; i++)
            aug.get(i).addAll(eye.get(i));

        var rref = rref(aug);

        // Check if the left half became I (otherwise singular)
        for (int i = 0; i < n; i++)
            if (Math.abs(rref.get(i).get(i) - 1.0) > EPS)
                throw new IllegalStateException("Matrix is singular and has no inverse.");

        // Extract the right half
        ArrayList<ArrayList<Double>> inv = new ArrayList<>();
        for (int i = 0; i < n; i++) {
            ArrayList<Double> row = new ArrayList<>();
            for (int j = n; j < 2 * n; j++)
                row.add(rref.get(i).get(j));
            inv.add(row);
        }
        return inv;
    }

    /**RANK**/

    /** Rank of a matrix = number of non-zero rows in its RREF. */
    public static int rank(ArrayList<ArrayList<Double>> matrix) {
        var r = rref(matrix);
        int count = 0;
        for (ArrayList<Double> row : r) {
            boolean nonZero = false;
            for (double v : row)
                if (Math.abs(v) > EPS) { nonZero = true; break; }
            if (nonZero) count++;
        }
        return count;
    }

    /**TRACE**/

    /** Trace = sum of main-diagonal elements. */
    public static double trace(ArrayList<ArrayList<Double>> A) {
        int n = Math.min(A.size(), A.get(0).size());
        double sum = 0;
        for (int i = 0; i < n; i++)
            sum += A.get(i).get(i);
        return sum;
    }

    /**PRINT**/

    /** Pretty-print a matrix with aligned columns. */
    public static void print(ArrayList<ArrayList<Double>> m) {
        for (ArrayList<Double> row : m) {
            for (double v : row)
                System.out.printf("%10.4f", v);
            System.out.println();
        }
        System.out.println();
    }
}


