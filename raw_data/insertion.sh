echo "intinate importing data $2  wiki at $1"

echo "Import catetory_link data"
mysql -u root --password=postech $2wiki < $2wiki-$1-categorylinks.sql
echo "Import page data"
mysql -u root --password=postech $2wiki < $2wiki-$1-page.sql
echo "Import lang-link data"
mysql -u root --password=postech $2wiki < $2wiki-$1-langlinks.sql
echo "$2wiki at $1 Cmpleted"
