# tutorial on youtube
# SQLite Databses With Python - Full Course
# freeCodeCamp.org

import sqlite3

# create new database
# run script in terminal to create the database
conn = sqlite3.connect('customer.db')

# to create a table, create a cursor first
c = conn.cursor()
# doc string """ """
c.execute("""CREATE TABLE customers (
		first_name text,
		last_name text,
		email text
	)""")

# 5 different datatypes:
# NULL exists?
# INTEGER number
# REAL decimals
# TEXT words
# BLOB images, music files

# put the stuff into the database
conn.commit()

# close connection
conn.close()




