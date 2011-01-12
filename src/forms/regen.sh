ls -d1 *.form | xargs -n1 ffc -l dolfin -f split_implementation
mv *.h ../unicorn
mv *.cpp ../
rm *.py
