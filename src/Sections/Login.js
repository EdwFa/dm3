import React, { Component } from 'react';

import { Navigate } from 'react-router-dom';
import { useState, useEffect } from 'react';

import { variables } from './Variables.js';


export class Login extends Component {

    constructor(props) {
        super(props);

        this.state = {
            token: variables.token,
            email: "",
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

      handleEmailChange = (event) => {
        this.setState({email: event.target.value});
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
          body: JSON.stringify({email: this.state.email, password: this.state.password}),
        })
        .then(this.isResponseOk)
        .then((data) => {
          console.log(data);
          this.setState({token: data.token.key});
          variables.token = data.token.key
          variables.email = data.email
        })
        .catch((err) => {
          console.log(err);
          this.setState({error: "Wrong email or password."});
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
                  <label htmlFor="email">Email</label>
                  <input type="email" className="form-control" id="email" name="email" value={this.state.email} onChange={this.handleEmailChange} />
                </div>
                <div className="form-group">
                  <label htmlFor="password">Password</label>
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
