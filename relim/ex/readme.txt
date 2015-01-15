The example files in this directory demonstrate how to use the
options -w, -r, -f, and -b. This file also explains the conversion
scripts, which can convert different input formats into the format
needed by the relim program.

In the file test1.tab transactions are separated by newline characters
and the items of a transaction are separated by spaces. This is the
standard input format and hence the file can be processed directly:
  relim test1.tab test1.out

In the file test2.tab the same transactions can be found, but several
different field separators are used. This file can be processed with:
  relim -f ",.;:" -l test2.tab test2.out

The file test3.tab has basically the same format as the file test1.tab,
with the only difference that the last fields of each record states an
(integer) transaction weight. This allows us to combine transactions,
so that test2.tab has only 8 lines, while test1.tab has 10 lines,
because the transactions "a b c" and "a b c d" occur twice. In order
to instruct the program to interpret the last field of each record as
such a weight, is has to be invoked with the option -w:
  relim -w test3.tab test3.out

The files test4.tab to test6.tab are in formats that may be common,
but which cannot be processed directly with the relim program.

In the file test4.tab each line contains a transaction identifier and
an item, separated by a space. This file can be converted into the
standard input format with the script tid2set, i.e., with
  tid2set test4.tab x.tab
Note that in this script the input file (here: test4.tab) is sorted
w.r.t. the transaction identifier, so that items belonging to the
same transaction occupy consecutive lines/records.

In the file test5.tab the first line states the item names and the
following lines contain flags T (true) and F (false) depending on
whether the item is contained in the transaction represented by the
line or not. This format can be converted into the standard input
format with the script flg2set, i.e., with
  flg2set test5.tab x.tab

In the file test5.tab there is one item per line and transactions
are separated by blank lines. This format can be converted into the
standard input format with the script row2set, i.e., with
  row2set test5.tab x.tab

The additional scripts tab2set and hdr2set convert tables with column
numbers or column names into a format appropriate for the relim
program. They are invoked in the same way as all other scripts
discussed above, i.e., with
  tab2set a.tab b.tab
or
  hdr2set a.tab b.tab
where a.tab is the name of the input file and b.tab the name of the
output file. The script tab2set replaces each table entry "x" of the
input file by "Xi=x", where i is the column number (starting with 1).
The script hdr2set reads the variable names from the first line of
the input file and then replaces each table entry "x" by "X=x", where
"X" is the variable name that was found in the corresponding column
of the first line. These scripts are handy if you want to process
tabular data by treating each table row as a transaction.

Note that any input may also be read from standard input and any output
may be sent to standard output, simply by specifying a '-' or an empty
string "" instead of a filename. For example
  cat test1.tab | relim - -
reads the transactions from standard input (where they are fed by the
cat command) and writes the item sets directly to the terminal. They
may be piped to any other program or script, since all other messages
of the relim program are written to standard error.

Enjoy,
Christian Borgelt
