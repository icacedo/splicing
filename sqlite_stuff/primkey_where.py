import sqlite3

# primary keys are a unique id number for each record

conn = sqlite3.connect('customer.db')

c = conn.cursor()
# select row id, and everything from customers
c.execute("SELECT rowid, * FROM customers")

items = c.fetchall()
for item in items:
	print(item)

print('#####')

# where clause
# can also use >=, >
#c.execute("SELECT * FROM customers WHERE last_name = 'Elder'")
# % is wildcard
c.execute("SELECT * FROM customers WHERE last_name LIKE 'Br%' ")

items = c.fetchall()

for item in items:
	print(item)


conn.commit()
conn.close()

