<?php
	//This function separates the extension from the rest of the file name and returns it
	$array = $_POST['algo'];
	if(empty($array)) die("You didn't select any consensus algorithm");

	//File automatically downloads as output.txt
	header('Content-Disposition: attachment; filename="output.txt"');

	//Build integer representing selected algorithms
	$N = count($array);
	$cnt = 0;
	for($i = 0; $i < $N; $i++) $cnt += (1 << $array[$i]);

	//The target file consists of the following parts:
	//1. The subdirectory you want to save into... make sure it exists!
	//2. A random number as the filename. You could also use a timestamp here if you prefer.
	//3. The extension, which we take to be .nex.
	$target = "upload/" . rand() . ".nex";

	// echo $_FILES['uploaded']['tmp_name'] . "<br />";
	if(move_uploaded_file($_FILES['uploaded']['tmp_name'], $target)){
		// echo "The file has been uploaded as " . $target . "<br>";
		$cmdLine = "echo " . $target . " " . $cnt . " " . $_POST["rooted"] . " | ./FACT";
		// echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);

		echo $output;

		$cmdLine = "rm " . $target;
		// echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
	}
	else echo "Sorry, there was a problem uploading your file.";
?>
