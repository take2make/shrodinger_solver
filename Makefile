NAME = shrodinger

CC = c++ -std=c++17 -I/usr/local/include/eigen3
CPP_FLAGS = -Wall -Wextra -Werror

OBJ_DIR = ./obj/
SRC_DIR = ./src/
INC_DIR = ./include/

SRC = $(wildcard $(SRC_DIR)*.cpp) $(wildcard $(SRC_DIR)*/*.cpp)
OBJ = $(SRC:$(SRC_DIR)%.cpp=$(OBJ_DIR)%.o)

all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(CPP_FLAGS) $(OBJ) -o $(NAME)

$(OBJ_DIR)%.o: $(SRC_DIR)%.cpp
	@mkdir -p $(dir $@)
	$(CC) $(CPP_FLAGS) -I$(INC_DIR) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR)

fclean: clean
	rm -r $(NAME)

re: fclean all

plot: 
	python3.11 scripts/plot_shrodinger.py

heatmap:
	python3.11 scripts/heatmap_proba.py

compare_1d:
	python3.11 scripts/check_gauss_1d.py

heatmap_2d:
	python3.11 scripts/heatmap_2D.py

compare_2d:
	python3.11 scripts/check_gauss_2d.py

proba:
	python3.11 scripts/check_proba.py

proba_2d:
	python3.11 scripts/check_proba_2D.py

test:
	./shrodinger && make heatmap

.PHONY: all clean fclean re