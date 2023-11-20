#!/bin/bash

echo "Start of the script"

sleep_duration=$((RANDOM % 10))

echo "Sleeping for $sleep_duration seconds..."
sleep $sleep_duration

echo "End of the script"