import java.util.Scanner;
public class hmm0 {
    public static void print_table_to_line(double[][] A){
        String out = String.format("%d %d ", A.length, A[0].length);
        for(int i=0;i<A.length;i++){
            for(int j=0;j<A[0].length;j++){
                out += String.format("%f ", A[i][j]);
            }
        }
        int tmp = out.length();
        System.out.println(out.substring(0, tmp-1));
    }
    public static double[][] multiplyTable(double[][] A, double[][] B){
        double[][] C = new double[A.length][B[0].length];
        for(int i=0;i<A.length;i++){
            for(int j=0;j<B[0].length;j++){
                double sum = 0;
                for (int k=0;k<A[i].length;k++){
                    sum+= A[i][k]*B[k][j];
                }
                C[i][j] = sum;
            }
        }
        return C;
    }
    public static void print_table(double[][] arr){
        for(int i=0;i<arr.length;i++){
            if (arr[i].getClass().isArray()){
                for(int j=0;j<arr[i].length;j++){
                    System.out.print(arr[i][j]);
                    System.out.print(' ');
                }
            }else{
                System.out.print(arr[i]);
            }
            System.out.println();
        }
    }
    public static double[][] readinputtable(String inp){
        Scanner s = new Scanner(inp);
        int rows = s.nextInt();
        int cols = s.nextInt();
        double[][] A = new double[rows][cols];
        for(int i=0;i<rows;i++){
            for(int j=0;j<cols;j++){
                A[i][j] = s.nextDouble();
            }
        }
        s.close();
        return A;
    }
    public static void main(String[] args) {
        Scanner scan = new Scanner(System.in);
        double[][] A = readinputtable(scan.nextLine());
        double[][] B = readinputtable(scan.nextLine());
        double[][] p = readinputtable(scan.nextLine());
        //print_table(A);
        //print_table(B);
        //print_table(p);
        double[][] tmp = multiplyTable(p,A);
        print_table_to_line(multiplyTable(tmp,B));
        scan.close();
    }
}
