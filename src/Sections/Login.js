import React, { Component } from 'react';

import { Navigate } from 'react-router-dom';
import { useState, useEffect } from 'react';

import { variables } from './Variables.js';


export class Login extends Component {

    constructor(props) {
        super(props);

        this.state = {
            token: variables.token,
            username: "",
            password: "",
            error: "",
            isAuthenticated: false,
        }
    }

    componentDidMount() {
        return
    }

    handlePasswordChange = (event) => {
        this.setState({password: event.target.value});
    }

      handleUserNameChange = (event) => {
        this.setState({username: event.target.value});
    }

    isResponseOk(response) {
        if (response.status >= 200 && response.status <= 299) {
          return response.json();
        } else {
          throw Error(response.statusText);
        }
    }

    login = (event) => {
        event.preventDefault();
        fetch(variables.API_URL + "/accounts/login", {
          method: "POST",
          headers: {
            "Content-Type": "application/json",
          },
          credentials: "same-origin",
          body: JSON.stringify({username: this.state.username, password: this.state.password}),
        })
        .then(this.isResponseOk)
        .then((data) => {
          console.log(data);
          this.setState({token: data.key});
          variables.token = data.key
        })
        .catch((err) => {
          console.log(err);
          this.setState({error: "Wrong username or password."});
        });
    }

    logout = () => {
        fetch(variables.API_URL + "/accounts/logout",
            {
              headers: {
                'Content-Type': 'application/json;charset=utf-8',
                'Authorization': `Token ${this.state.token}`,
              },
            }
        )
        .then(this.isResponseOk)
        .then((data) => {
          console.log(data);
          this.setState({token: null});
          variables.token = null;
        })
        .catch((err) => {
          console.log(err);
        });
    };

    render() {
        if (!this.state.token) {
          return (
            <div className="container mt-3">
              <h1>React Cookie Auth</h1>
              <br />
              <h2>Login</h2>
              <form onSubmit={this.login}>
                <div className="form-group">
                  <label htmlFor="username">Username</label>
                  <input type="text" className="form-control" id="username" name="username" value={this.state.username} onChange={this.handleUserNameChange} />
                </div>
                <div className="form-group">
                  <label htmlFor="username">Password</label>
                  <input type="password" className="form-control" id="password" name="password" value={this.state.password} onChange={this.handlePasswordChange} />
                  <div>
                    {this.state.error &&
                      <small className="text-danger">
                        {this.state.error}
                      </small>
                    }
                  </div>
                </div>
                <button type="submit" className="btn btn-primary">Login</button>
              </form>
            </div>
          );
        } else {
          return <Navigate push to="/" />
        }
    }
}
