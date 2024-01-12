import sqlite3

# make functions to be called from the other file
# query the db and return all records
def show_all():
	
	conn = sqlite3.connect('customer.db')

	c = conn.cursor()

	c.execute("SELECT rowid, * FROM customers")
	items = c.fetchall()
	for item in items:
		print(item)

	conn.commit()
	conn.close()

# add a new record to the table
def add_one(first, last, email):
	
	conn = sqlite3.connect('customer.db')
	c = conn.cursor()
	c.execute("INSERT INTO customers VALUES (?,?,?)", (first, last, email))	
	conn.commit()
	conn.close()

# delete record from the table
def delete_one(id):
	conn = sqlite3.connect('customer.db')
	c = conn.cursor()
	c.execute("DELETE from customers WHERE rowid = (?)", id)

	conn.commit()
	conn.close()

# create a function to add many records
def add_many(list):
	conn = sqlite3.connect('customer.db')
	c = conn.cursor()
	c.executemany("INSERT INTO customers VALUES (?,?,?)", (list))

	conn.commit()
	conn.close()

# lookup with where
def email_lookup(email):
	conn = sqlite3.connect('customer.db')
	c = conn.cursor()
	c.execute("SELECT * from customers WHERE email = (?)", (email,))

	items = c.fetchall()
	for item in items:
		print(item)

	conn.commit()
	conn.close()








