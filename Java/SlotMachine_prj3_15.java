import java.util.*;
//by Alex Fitz

public class SlotMachine_prj3_15
{
   public static void main (String[] arg)
   {
      Random generator = new Random();
      int num1, num2, num3;
      String yes="y";
      Scanner scan = new Scanner (System.in);
      
      while (yes.equalsIgnoreCase("y")) //allows y or Y
      {     
         num1 = generator.nextInt(10);
         System.out.print (num1 + " ");//prints first random number
      
         num2 = generator.nextInt(10);
         System.out.print (num2 + " ");//prints second random number
      
         num3 = generator.nextInt(10);
         System.out.print (num3 + " ");//prints third random number
         
         System.out.println ();
         if (num1 == num2 && num2 == num3)
            System.out.println ("All three numbers are the same");
         else if (num1 == num2||num2==num3||num1==num3)
            System.out.println ("two numbers are the same");
         else
            System.out.println("The numbers are different from each other");
      
              
         System.out.println();
      
         if (yes.equalsIgnoreCase("y"))
         {
            System.out.print ("Play again (Y to play again)? ");
            yes = scan.nextLine();
         }
      
         
      }
      System.out.println ("Play again soon!");
   
   
   
   }
}