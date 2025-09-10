# Create a README file
echo "# particle_aggregation" >> README.md

# Initialize a new Git repository
git init

# Add the README file to the staging area
git add README.md

# Commit the changes
git commit -m "first commit"

# Create and switch to the main branch
git branch -M main

# Add your remote repository (change the URL to your repository)
git remote add origin https://github.com/Turjokez/particle_aggregation.git

# Push the changes to GitHub
git push -u origin main
