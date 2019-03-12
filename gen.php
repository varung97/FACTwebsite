<?php require_once "top.php"; ?>

<form enctype="multipart/form-data" align = 'left' action="upload_file.php" method="POST">
Please choose a file:<br> <input name="uploaded" type="file" /><br /><br><br>
Please select the type of consensus algorithm you want to use. Note: For complexity analysis, k = number of trees, n = number of leaves<br><br><br>
Click <a href = 'sample.nex'>here</a> for a sample input file.
<br><br>
<input type="checkbox" name="algo[]" value="0">Strict Consensus: O(kn)<br>
<br>
<input type="checkbox" name="algo[]" value="1">Slow Majority Consensus: O(k<sup>2</sup>n)<br>
<input type="checkbox" name="algo[]" value="2">Fast Majority Consensus V1: O(kn lg k)<br>
<input type="checkbox" name="algo[]" value="3">Fast Majority Consensus V2: O(kn)<br>
<br>
<input type="checkbox" name="algo[]" value="4">Slow Greedy Consensus:   O(kn<sup>2</sup> + kn<sup>3</sup> + n<sup>2</sup>)<br>
<input type="checkbox" name="algo[]" value="5">Fast Greedy Consensus:   O(kn<sup>2</sup>)<br>
<br>
<input type="checkbox" name="algo[]" value="6">Slow Loose Consensus: O(k<sup>2</sup>n<sup>2</sup>)<br>
<input type="checkbox" name="algo[]" value="7">Fast Loose Consensus: O(kn)<br>
<br>
<input type="checkbox" name="algo[]" value="8">Majority Rule(+) Consensus: O(kn)<br>
<input type="checkbox" name="algo[]" value="9">Adams Consensus: O(kn<sup>2</sup>)<br>
<input type="checkbox" name="algo[]" value="10">Adams Consensus Fast: O(kn lg n)<br>
<br>


Please select whether the input trees are rooted:<br>
<input type="radio" name="rooted" value="1" checked>Rooted<br>
<input type="radio" name="rooted" value="0">Unrooted<br><br><br>
<input type="submit" value="Upload" />
</form> 
<?php require_once "bottom.php"; ?>
