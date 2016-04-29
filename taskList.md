#### Task list, week of April 29-May 6

For database queries, use the Azure spatial_database. Below I use __bold__ text to represent table names and _italic_ text to represent variables in the table or variables in the Python script. 

- [ ] Look at variable _masconGeoms_. It is a Pandas dataframe generated from a SQL query to the database. It contains the _mascon_ number and the WKT geometry of each square polygon. I've already used the _geom_ field to generate a shapely object called _masconShapely_.TASK: add label of _mascon_ number to the cartopy plot. Use the centroid method on the _masconShapely_ object.   
- [ ] I need to know how much glacier ice is in each _mascon_ geometry. TASK: write a SQL query to generate a sum of _area_ from __modern__ GROUP BY _geom_ of __mascon_fit__. Write this query so that it returns a dataframe with two columns: _mascon_, _sumGlacArea_.