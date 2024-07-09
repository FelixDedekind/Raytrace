import cv2
import os
import pandas as pd

# Initialize global variables
points = []
distances = []

# Mouse callback function to capture points
def capture_points(event, x, y, flags, param):
    global points
    global distances
    if event == cv2.EVENT_LBUTTONDOWN:
        points.append((x, y))
        cv2.circle(img, (x, y), 5, (0, 255, 0), -1)
        cv2.imshow('Image', img)
        if len(points) == 2:
            distance = int(((points[1][0] - points[0][0]) ** 2 + (points[1][1] - points[0][1]) ** 2) ** 0.5)
            print(f"distance is: {distance}")
            distances.append(distance)
            points = []

# Load an image

# Define the directory path
directory_path = 'images'
output_path = 'magnification/out/out1.csv'

# Loop through each file in the directory

for filename in os.listdir(directory_path):
    file_path = os.path.join(directory_path, filename)
    if os.path.isfile(file_path):  # Check if it's a file
        img_path = 'images/65.jpg'
        img = cv2.imread(img_path)
        cv2.imshow('Image', img)
        cv2.setMouseCallback('Image', capture_points)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

data = {
    'distances': distances,
}
df = pd.DataFrame(data)

df.to_csv(output_path)
        

    