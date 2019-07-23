#ifndef SHADER_H
#define SHADER_H

#include <GL/glew.h>
#include <glm/glm.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Shader
{
public:
	Shader(const char* vsPath, const char* fsPath, const char* gsPath=NULL)
	{
		std::string vertexShaderCode;
		std::string fragmentShaderCode;
		std::string geometryShaderCode;

		std::ifstream vertexShaderFile;
		std::ifstream fragmentShaderFile;
		std::ifstream geometryShaderFile;

		vertexShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		fragmentShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);
		geometryShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

		try 
		{
			vertexShaderFile.open(vsPath);
			fragmentShaderFile.open(fsPath);
			std::stringstream vertexShaderStream, fragmentShaderStream;
			
			vertexShaderStream << vertexShaderFile.rdbuf();
			fragmentShaderStream << fragmentShaderFile.rdbuf();
			
			vertexShaderFile.close();
			fragmentShaderFile.close();

			vertexShaderCode = vertexShaderStream.str();
			fragmentShaderCode = fragmentShaderStream.str();
			
			if (gsPath != NULL)
			{
				geometryShaderFile.open(gsPath);
				std::stringstream geometryShaderStream;
				geometryShaderStream << geometryShaderFile.rdbuf();
				geometryShaderFile.close();
				geometryShaderCode = geometryShaderStream.str();
			}
		}
		catch (std::ifstream::failure error)
		{
			std::cout << "SHADER READING FAILED! " << std::endl;
		}

		const char* vsCode = vertexShaderCode.c_str();
		const char* fsCode = fragmentShaderCode.c_str();

		unsigned int vs, fs, gs;
		vs = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(vs, 1, &vsCode, NULL);
		glCompileShader(vs);
		checkCompileError(vs, "VERTEX");

		fs = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(fs, 1, &fsCode, NULL);
		glCompileShader(fs);
		checkCompileError(fs, "FRAGMENT");

		if (gsPath != NULL)
		{
			const char* gsCode = geometryShaderCode.c_str();
			gs = glCreateShader(GL_GEOMETRY_SHADER);
			glShaderSource(gs, 1, &gsCode, NULL);
			glCompileShader(gs);
			checkCompileError(gs, "GEOMETRY");
		}

		ID = glCreateProgram();
		glAttachShader(ID, vs);
		glAttachShader(ID, fs);
		if (gsPath != NULL) glAttachShader(ID, gs);
		glLinkProgram(ID);
		checkCompileError(ID, "PROGRAM");

		glDeleteShader(vs);
		glDeleteShader(fs);
		if (gsPath != NULL) glDeleteShader(gs);
	}

	void use()
	{
		glUseProgram(ID);
	}

	// utility uniform functions
	// ------------------------------------------------------------------------
	void setBool(const std::string &name, bool value) const
	{
		glUniform1i(glGetUniformLocation(ID, name.c_str()), (int)value);
	}
	// ------------------------------------------------------------------------
	void setInt(const std::string &name, int value) const
	{
		glUniform1i(glGetUniformLocation(ID, name.c_str()), value);
	}
	// ------------------------------------------------------------------------
	void setFloat(const std::string &name, float value) const
	{
		glUniform1f(glGetUniformLocation(ID, name.c_str()), value);
	}
	// ------------------------------------------------------------------------
	void setVec2(const std::string &name, const glm::vec2 &value) const
	{
		glUniform2fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
	}
	void setVec2(const std::string &name, float x, float y) const
	{
		glUniform2f(glGetUniformLocation(ID, name.c_str()), x, y);
	}
	// ------------------------------------------------------------------------
	void setVec3(const std::string &name, const glm::vec3 &value) const
	{
		glUniform3fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
	}
	void setVec3(const std::string &name, float x, float y, float z) const
	{
		glUniform3f(glGetUniformLocation(ID, name.c_str()), x, y, z);
	}
	// ------------------------------------------------------------------------
	void setVec4(const std::string &name, const glm::vec4 &value) const
	{
		glUniform4fv(glGetUniformLocation(ID, name.c_str()), 1, &value[0]);
	}
	void setVec4(const std::string &name, float x, float y, float z, float w)
	{
		glUniform4f(glGetUniformLocation(ID, name.c_str()), x, y, z, w);
	}
	// ------------------------------------------------------------------------
	void setMat2(const std::string &name, const glm::mat2 &mat) const
	{
		glUniformMatrix2fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}
	// ------------------------------------------------------------------------
	void setMat3(const std::string &name, const glm::mat3 &mat) const
	{
		glUniformMatrix3fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}
	// ------------------------------------------------------------------------
	void setMat4(const std::string &name, const glm::mat4 &mat) const
	{
		glUniformMatrix4fv(glGetUniformLocation(ID, name.c_str()), 1, GL_FALSE, &mat[0][0]);
	}

private:
	void checkCompileError(GLuint shader, std::string shaderType)
	{
		GLint success;
		GLchar infoLog[1024];

		if (shaderType != "PROGRAM")
		{
			glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
			if (!success)
			{
				glGetShaderInfoLog(shader, 1024, NULL, infoLog);
				std::cout << shaderType << " SHADER COMPILE FAILED: " << infoLog << std::endl;
			}
		}
		else
		{
			glGetProgramiv(shader, GL_LINK_STATUS, &success);
			if (!success)
			{
				glGetProgramInfoLog(shader, 1024, NULL, infoLog);
				std::cout << shaderType << " LINK FAILED: " << infoLog << std::endl;
			}
		}
	}

public:
	unsigned int ID;
};


#endif