<?php include_once "top.php"; ?>
<p>This is a guide!</p>

<p>The program accepts files from 10ktrees.</p>

<p>The program also accepts files of the following format: </p>

<pre>
#NEXUS
BEGIN TREES;
1 A,
2 B,
3 C,
4 D,
5 E,
6 F,
7 G,
8 H,
9 I,
10 J,
11 K;
tree 1 = (8,(((9,4),((1,5),((7,6),3))),(2,(10,11))));
tree 2 = (8,((((1,5),((6,7),3)),(4,9)),(2,(11,10))));
tree 3 = (8,((2,(11,10)),((9,4),((1,5),(6,(7,3))))));
END;



</pre>

<p>Note: If the tree is unrooted, then the tree must not have any node of degree 2.</p>
<?php include_once "bottom.php"; ?>