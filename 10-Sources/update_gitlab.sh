echo "Gitlab Start"

cd Ds2D0enue_Analysis
# git switch KEK_Server
git status
git add .
git commit -m "update"
# git config --global user.name "adilmub"
# git config --global user.email amubarak@iastate.edu
echo "Ready to Push" 
git push -uf origin main

# Username: adilmub, Password: AzurBLUE8978
# git config credential.helper store 

cd ~
echo "Push Succesful"

rm /home/belle2/amubarak/Ds2D0enue_Analysis/10-Sources/update_gitlab.sh -f
cp update_gitlab.sh /home/belle2/amubarak/Ds2D0enue_Analysis/10-Sources
echo "Updated Source File"