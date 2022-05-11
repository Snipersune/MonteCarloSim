#!/bin/bash

ls data/008_* > data/008all_files
ls data/016_* > data/016all_files
ls data/032_* > data/032all_files
ls data/064_* > data/064all_files

VAR8=$(echo $(cat data/008all_files))
VAR16=$(echo $(cat data/016all_files))
VAR32=$(echo $(cat data/032all_files))
VAR64=$(echo $(cat data/064all_files))

tail -q -n 1 ${VAR8} > data/_008_data
tail -q -n 1 ${VAR16} > data/_016_data
tail -q -n 1 ${VAR32} > data/_032_data
tail -q -n 1 ${VAR64} > data/_064_data

rm data/008all_files
rm data/016all_files
rm data/032all_files
rm data/064all_files