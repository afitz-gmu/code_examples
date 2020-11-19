import java.util.*;


public class RockPaperScissors
{
   public static void main (String[] args)
   {
   
      String userMove = ".";
      String cpuMove = "";
      
      int cpuNum;
      int lose = 0;
      int win = 0;
      int tie = 0;
      int userNum;
     
   
      Scanner scan = new Scanner(System.in);
      Random generator = new Random();
   
      System.out.println ("Let's play Rock, Paper, Scissors! Rock = R, Paper = P, and Scissors = S.");
      System.out.print ("Enter your move. Enter stop to quit: ");
         
      while (!userMove.equals("STOP"))
        
      { 
         userMove = scan.nextLine();
         userMove = userMove.toUpperCase();  //to make the comparison between the cpu generated number and the user input number
      
         cpuNum = generator.nextInt(3)+1;
         if   (cpuNum == 1)
            cpuMove = "R";
         else if (cpuNum == 2) 
            cpuMove = "S";
         else if (cpuNum == 3 )
            cpuMove = "P";
           
         
         
         if (cpuMove.equals(userMove))
         {
            System.out.println ("It's a tie");
            tie++;
         }
         else if (userMove.equals("R"))
         {
            if (cpuMove.equals("S"))
            {
               System.out.println ("The rock beats the scissors, you win");
               win++;
            }
            else if (cpuMove.equals("P"))
            {
               System.out.println ("The paper beats the rock, you lose");
               lose++;
            }
         }
         else if (userMove.equals("S"))
         {
            if (cpuMove.equals("R"))
            {
               System.out.println ("The rock beats the scissors, you lose");
               lose++;
            }
            else if (cpuMove.equals("P"))
            {
               System.out.println ("The scissors beat the paper, you win");
               win++;
            }
         }
         else if (userMove.equals("P"))
         {
            if (cpuMove.equals("R"))
            {
               System.out.println ("The paper beats the rock, you win");
               win++;
            }
            else if (cpuMove.equals("S"))
            {
               System.out.println("The Scissors beat the paper, you lose");
               lose++;
            }
         }
         else System.out.println ("Your move was invalid");
         
                           
         System.out.println ("The computer chose: " + cpuMove);
         System.out.println ("You chose: " +userMove);
         
         System.out.print ("Guess again: ");
         
      }
      
      {
         System.out.println ("You quit");
         System.out.print ("You had " + win + " wins, " + tie + " ties, and " + lose + " loses ");
      }
       
      
   }
}
   
