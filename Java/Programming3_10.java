import java.util.Scanner;
import java.util.Random;

public class Programming3_10
{
   public static void main (String[] args)
   {
      int answer, guess;
      int count = 1;
      int input = 0;
      final int HIGH = 100;
      final int SENTINEL = 999;
      Scanner scan = new Scanner (System.in);
      String yes="y";
   
      Random generator = new Random ();
      answer = generator.nextInt(HIGH)+1;
   
            
      System.out.print ("I'm thinking of a number between 1 through 100. Guess what that number is, and if you would like to quit enter 999: ");
      guess = scan.nextInt();
      
   
      while (yes.equalsIgnoreCase("y"))
      {
         while (guess != answer && guess != SENTINEL)
         {
            count++;
                  
            if (guess > answer)
               System.out.print ("Your guess was too high, guess again: ");               
               
            else 
               System.out.print ("Your guess was too low, guess again: ");
            guess = scan.nextInt();         
         }
      
      
         if (guess == SENTINEL)
         {
            System.out.println ("You quit");
            System.out.print ("Times it took you to guess the correct number: " + count);
         }
         else
         {
            System.out.println ("YOU GUESSED CORRECTLY!");
            System.out.println ("Times it took you to guess the correct number: " + count);
         }
         
         System.out.println();
         
         
         System.out.println ("Play again (Y to play again)? ");
         yes = scan.nextLine();
         
      }
   }
}
