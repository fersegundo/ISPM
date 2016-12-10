<?php
 
// sleep for 2 sec show that the androd swipe refresh will be visible for sometime
//sleep(2);
 
// "DATABASE"
$all_names = array(             "Let*s bi Chicks!",
                    "National Geographic Channelllllllllx",
                    "Dj* Workshop",
                    "Klausurtag Zukunftswerkstatt und Jahresplanung 2017",
                    "Discovery",
                    "Ramdi Tamda",
                    "KTM",
                    "Harley Davison",
                    "Rock and Girl",
                    "Bla bla bla",
                    "Otra movida",
                    "Algo mas",
                    );
$all_images = array(            "https://www.lambda-bb.de/wp-content/uploads/2012/12/20150218_195720-300x225.jpg",
                    "https://www.lambda-bb.de/wp-content/uploads/2012/12/DSC04866-300x225.jpg",
                    "https://www.lambda-bb.de/wp-content/uploads/2012/12/IMG_20160817_192933-300x168.jpg",
                    null,
                    "http://api.androidhive.info/feed/img/discovery_mos.jpg",
                    "http://api.androidhive.info/feed/img/nav_drawer.jpg",
                    "http://api.androidhive.info/feed/img/ktm_1290.jpg",
                    "http://api.androidhive.info/feed/img/harley_bike.jpg",
                    "http://api.androidhive.info/feed/img/rock.jpg",
                    null,
                    "http://api.androidhive.info/feed/img/imagenquenoexiste.jpg",
                    "http://cdn.slashgear.com/wp-content/uploads/2012/07/asdfdsa1.png",
                    );
$all_status = array( "Wir sind die LeT*s Bi Chicks, was fur lesbisch, bi und trans* steht. Vor allem sind wir aber eine Gruppe queerer, bunter Weiblichkeit_en. Wir sind alle zwischen 14 und 27 Jahren alt und kommen aus unterschiedlichsten Kontexten: Manche von uns arbeiten, manche studieren, andere gehen zur Schule, machen Praktika oder schauen noch.",
                    "Was genau machen diese Dj*s eigentlich wenn sie am Pult stehen? Welche Knopfe drehen, schieben, drucken sie wann und warum und wie machen sie das mit diesen ubergangen? Dieser Workshop bietet eine Einfuhrung zu dem Thema Auflegen mit MP3. Ich erklare euch grundlegend die Technik mit der ihr es zu tun haben werdet. Ausserdem gibt es eine Einfuhrung in das Programm Traktor. Bitte bringt die Musik mit, mit der ihr umprobieren wollt. (10-15 Songs als Mp3)",
                    "Wir planen das Jahr 2017. Bringt eure Ideen mit.",
			"muerto",
                    "disponible",
                    "no esta en su casa",
                    "que es esto?",
                    "quibra",
                    "ver imagen",
                    "no se que poner aqui",
                    "imagen erronea",
                    "imagen de internet cualquiera",
                    );
$all_profile_pic = array( "https://media.licdn.com/mpr/mpr/shrinknp_200_200/AAEAAQAAAAAAAANkAAAAJDI0ZDJiN2E3LTVmNmMtNGM4Ni1hZTk4LTgyODFkNjgwNTBmNg.jpg",
                    "http://api.androidhive.info/feed/img/nat.jpg",
                    "http://api.androidhive.info/feed/img/time.png",
                    "http://api.androidhive.info/feed/img/lincoln.jpg",
                    "http://api.androidhive.info/feed/img/discovery.jpg",
                    "http://api.androidhive.info/feed/img/ravi_tamada.jpg",
                    "http://api.androidhive.info/feed/img/ktm.png",
                    "http://api.androidhive.info/feed/img/harley.jpg",
                    "http://api.androidhive.info/feed/img/rock_girl.jpg",
                    "http://api.androidhive.info/feed/img/gandhi.jpg",
                    null,
                    "https://upload.wikimedia.org/wikipedia/commons/4/41/Cicero.PNG",
                    );
$all_timestamp = array( "1403375851950",
                    "1403375851950",
                    "1403375851950",
                    "1403375851950",
                    "1403375851960",
                    "1403375851970",
                    "1403375851980",
                    "1403375851990",
                    "1403375852000",
                    "1403375852100",
                    "1403375852200",
                    "1403375852300",
                   );
$all_url = array( null,
                    null,
                    "http://ti.me/1qW8MLB",
                    null,
                    "http://dsc.tv/xmMxD",
                    "http://www.androidhive.info/2013/11/android-sliding-menu-using-navigation-drawer/",
                    "",
                    "http://bit.ly/1wmBWaN",
                    "www.google.com",
                    "http://estonoexistenideconia.no/nideconia",
                    "y tambien se puede poner cualquier cosa",
                    "cualquierweb.com",
                   );
// reading offset from get parameter
$offset = isset($_GET['offset']) && $_GET['offset'] != '' ? $_GET['offset'] : 0;

// id of newest item
$newest_id = 12;

// page limit
$limit = 5;
 
 
$feed_array = array();

// loop through all
for ($j = $offset; $j < $offset + $limit && $j < sizeof($all_names); $j++) {
    $tmp = array();
    $tmp['id'] = $newest_id - $j;
    $tmp['name'] = $all_names[$j];
    $tmp['image'] = $all_images[$j];
    $tmp['status'] = $all_status[$j];
    $tmp['profilePic'] = $all_profile_pic[$j];
    $tmp['timeStamp'] = $all_timestamp[$j];
    $tmp['url'] = $all_url[$j];

    array_push($feed_array, $tmp);
}
$feed_output = array();
$feed_output['feed'] = $feed_array;
// printing json response
$json = json_encode($feed_output,  JSON_PRETTY_PRINT | JSON_UNESCAPED_SLASHES );
header('Content-Type: application/json');
echo $json;
//printf("%s", $json);
//printf("<pre>%s</pre>", $json);

?>