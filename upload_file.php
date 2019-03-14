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

	//We use the extension nex for our input file
	$ext = 'nex';

	//This line assigns a random number to a variable. You could also use a timestamp here if you prefer.
	$ran = rand();

	//This takes the random number (or timestamp) you generated and adds a . on the end, so it is ready of the file extension to be appended.
	$ran2 = $ran . ".";

	//This assigns the subdirectory you want to save into... make sure it exists!
	$dir = "upload/";
	//This combines the directory, the random file name, and the extension
	$target = $dir . $ran2 . $ext;
	$inputParams = $dir . $ran2 . "in";

	// echo $_FILES['uploaded']['tmp_name'] . "<br />";
	if(move_uploaded_file($_FILES['uploaded']['tmp_name'], $target)){
		// echo "The file has been uploaded as " . $target . "<br>";
		$cmdLine = "echo " . $target . " " . $cnt . " " . $_POST["rooted"] . " > " . $inputParams;
		// echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine,"r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
		// echo "finished executing: " . $cmdLine . "<br />";
		$cmdLine = "./FACT < " . $inputParams;
		// echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);

		echo $output;

		$cmdLine = "rm " . $inputParams . "; rm " . $target;
		// echo "executing: " . $cmdLine . "<br />";
		$outhandle = popen($cmdLine, "r");
		$output = stream_get_contents($outhandle,-1,-1);
		pclose($outhandle);
	}
	else echo "Sorry, there was a problem uploading your file.";
?>
