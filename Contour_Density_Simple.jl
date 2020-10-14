using Plotly

size = 100
x = range(-2*pi, stop=2*pi, length=size)
y = range(-2*pi, stop=2*pi, length=size)
z = rand(size, size)
for i = 1:size
  for j = 1:size
    r2 = (x(i)^2 + y(j)^2)
    z(i,j) = sin(x(i))*cos(y(j))*sin(r2)/log(r2+1)
  end
end

data = [
  [
    "z" => z,
    "x" => x,
    "y" => y,
    "type" => "contour"
  ]
]
response = Plotly.plot(data, ["filename" => "simple-contour", "fileopt" => "overwrite"])
plot_url = response["url"]
