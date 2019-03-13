<?php
	//This function separates the extension from the rest of the file name and returns it
	$array = $_POST['algo'];
	if(empty($array)) die("You didn't select any consensus algorithm");
	function findexts ($filename) {
		$filename = strtolower($filename) ;
		$exts = split("[/\\.]", $filename) ;
		$n = count($exts)-1;
		$exts = $exts[$n];
		return $exts;
	}
	$N = count($array);
	$cnt = 0;
	for($i = 0;$i < $N;$i++) $cnt += (1<<$array[$i]);
	// TODO: What's the point? Why not just use .nex?
	//This applies the function to our file
	$ext = findexts ($_FILES['uploaded']['name']) . 'nex';

	//This line assigns a random number to a variable. You could also use a timestamp here if you prefer.
	$ran = rand () ;

	//This takes the random number (or timestamp) you generated and adds a . on the end, so it is ready of the file extension to be appended.
	$ran2 = $ran.".";

	//This assigns the subdirectory you want to save into... make sure it exists!
	$dir = "upload/";
	//This combines the directory, the random file name, and the extension
	$target = $dir . $ran2 . $ext;

	echo $target . "<br />";

	echo $_FILES['uploaded']['tmp_name'] . "<br />";
	if(move_uploaded_file($_FILES['uploaded']['tmp_name'], $target)){
		echo "The file has been uploaded as ".$ran2.$ext."<br>";
		$cmdLine = "echo ".$target." ".$cnt." ".$_POST["rooted"]." > ".$dir.$ran2."in";
		echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine,"r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
		echo "finished executing: " . $cmdLine . "<br />";
		$cmdLine = "./a < " . $dir.$ran2."in" . " > " . $dir . $ran2 . "txt";
		//$cmdLine = "a < " . $dir.$ran2."in" . " > " . $dir . $ran2 . "txt";
		echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
		echo "finished executing: " . $cmdLine . "<br />";
		$cmdLine = "chmod 755 " . $dir . $ran2 . "txt";
		echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
		echo "<a href = '".$dir.$ran2."txt"."'>Download Output</a><br>";
	}
	else echo "Sorry, there was a problem uploading your file.";
?>
