#!/bin/bash
echo "Hello World !"
name='xuheng'
echo $name
echo "$name"
echo \{$name\}
for skill in Ada Coffe Action Jave;do
    echo "I am good at ${skill}Script"
done

your_name="qinjx"
greeting="hello, "$your_name" !"
greeting_1="hello, ${your_name} !"
echo $greeting $greeting_1
