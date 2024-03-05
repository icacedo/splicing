import sqlite3

conn = sqlite3.connect('customer.db')

c = conn.cursor()

# update records
# this would change all records where last name is Elder
# use rowid to be specific
c.execute("""UPDATE customers SET first_name = 'Bob'
			WHERE rowid = 1
	""")

conn.commit()

# delete or drop records
#c.execute("DELETE from customers WHERE rowid = 6")

# order by 
# query database
#c.execute("SELECT rowid, * FROM customers ORDER BY rowid DESC")
#c.execute("SELECT rowid, * FROM customers ORDER BY last_name")

# and/or

#c.execute("SELECT rowid, * FROM customers WHERE last_name Like 'Br%' AND rowid = 3")
# can have many OR statements
#c.execute("SELECT rowid, * FROM customers WHERE last_name Like 'Br%' OR rowid = 3")

# limit results
#c.execute("SELECT rowid, * FROM customers LIMIT 2")
#c.execute("SELECT rowid, * FROM customers ORDER BY rowid DESC Limit 2")

# delete table
#c.execute("DROP TABLE customers")

items = c.fetchall()

for item in items:
	print(item)

conn.close()










