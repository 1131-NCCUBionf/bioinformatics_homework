# pro2. calculate the Sum-Of-Pair score of the MSA

-   your name
-   student ID

## Description

-   Write a Python program to calculate the multiple sequence alignment's sum-of-pair score (SoP).
-   Creating your own program, i.e., hw2.py.
-   Packages you can use: numpy, pandas
-   You write a program with a function named calculate_SoP, ie.

```
def calculate_SoP(input_path, score_path, gopen, gextend):
    .
    .
    .
    .
```

## File

-   hw2_ref.py: You can start from this reference code and try to write your own comment in English
-   pam100.txt
-   pam250.txt
-   test1.fasta

## Parameters

-   input_path: fasta file (ex. test1.fasta)
-   score_path: score file (ex. pam250.txt)
-   gopen: gap open penalty
-   gextend: gap extend penalty

## Command

Please go ahead and execute your code using the following command.

```Python
calculate_SoP("examples/test1.fasta", "pam250.txt", -10, -2) #score=1047
calculate_SoP("examples/test2.fasta", "pam100.txt", -8, -2) #score=606
```

## Evaluation

10 testing data(5 public, 5 private)

The correct answer gets 10 points for each testing data.

### Penalty

-   High code similarity to others: YOUR SCORE = 0

## References

https://chatgpt.com/share/67259d48-49ec-800a-8d5a-7f2e3286ed47
