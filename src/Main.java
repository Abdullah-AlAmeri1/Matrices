
import java.util.*;



public class Main {
    public static void main(String[] args) {

        // 3×3 test matrix
        ArrayList<ArrayList<Double>> A = new ArrayList<>(List.of(
                new ArrayList<>(List.of(2.0, 1.0, -1.0)),
                new ArrayList<>(List.of(-3.0, -1.0, 2.0)),
                new ArrayList<>(List.of(-2.0, 1.0, 2.0))
        ));

        System.out.println("=== Original Matrix A ===");
        Matrix.print(A);

        System.out.println("=== RREF ===");
        Matrix.print(Matrix.rref(A));

        System.out.println("=== REF (upper triangular) ===");
        Matrix.print(Matrix.ref(A));

        System.out.println("=== Transpose ===");
        Matrix.print(Matrix.transpose(A));

        System.out.println("=== Inverse ===");
        Matrix.print(Matrix.inverse(A));

        System.out.println("=== Determinant ===");
        System.out.printf("det(A) = %.4f%n%n", Matrix.determinant(A));

        System.out.println("=== Rank ===");
        System.out.println("rank(A) = " + Matrix.rank(A) + "\n");

        System.out.println("=== Trace ===");
        System.out.printf("trace(A) = %.4f%n%n", Matrix.trace(A));

        System.out.println("=== LU Decomposition ===");
        var lu = Matrix.decomposeLU(A);
        System.out.println("L:"); Matrix.print(lu.L);
        System.out.println("U:"); Matrix.print(lu.U);

        System.out.println("=== LDU Decomposition ===");
        var ldu = Matrix.decomposeLDU(A);
        System.out.println("L:"); Matrix.print(ldu.L);
        System.out.println("D:"); Matrix.print(ldu.D);
        System.out.println("U:"); Matrix.print(ldu.U);

        System.out.println("=== Solve A·x = b ===");
        ArrayList<Double> b = new ArrayList<>(List.of(8.0, -11.0, -3.0));
        var sol = Matrix.solve(A, b);
        sol.forEach((k, v) -> System.out.printf("%s = %.4f%n", k, v));
    }
}

