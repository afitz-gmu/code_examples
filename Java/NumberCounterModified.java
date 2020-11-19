import java.util.Scanner;

//by Alex Fitz

public class NumberCounterModified
{
   public static void main (String [] args)
   {
   
      final int RANGE = 51;
      final int N_VALUES = 10;
      Scanner input = new Scanner(System.in);
      int counter = 1;
      
      int[] entered = new int[RANGE];
      
            
      for (int count = 0; count < N_VALUES; count++)
      
      {
      
         System.out.print ("Enter number " + counter + " (-25 to 25 inclusive): ");
         int value = input.nextInt();
         entered[value+25]++;  
         counter++;
         
      
      }
      
      {
         System.out.println ();
         
         for (int number = 0; number < RANGE; number++)
         {
            System.out.println((number-25) + ": " +entered[number]);
         }
        
      }
   }
}
